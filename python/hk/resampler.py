#!/usr/bin/env python3

import bisect
import so3g
from spt3g import core
import numpy as np
from scipy import interpolate

# Unused so far.
_G3VECTOR_NUMERIC_TYPES = [
    core.G3VectorComplexDouble,
    core.G3VectorInt,
    core.G3VectorTime,
]

def _scipy_interpolate(x_in, y_in, x_out, **kwargs):
    # Use the scpiy interpolator, with default settings.

    # Numpy-ize the x-axis to which data are to be interpolated.
    xx_out = np.array([t.time for t in x_out])
    if xx_out.shape[0] == 0:
        return np.array([])

    # Figure out which portion of the input vector overlaps with the output
    # range, as well as the index just beyond the output range.
    in_range = np.where((x_in[0] <= xx_out) & (x_in[-1] >= xx_out))[0]
    next_x_in = bisect.bisect_left(x_in, xx_out[-1])
    if next_x_in >= x_in.shape[0]:
        next_x_in = x_in.shape[0] - 1

    if len(in_range) > 1:
        # Should be the ordinary case: at least two points to interpolate.
        f = interpolate.interp1d(x_in, y_in, assume_sorted=True, **kwargs)
        front = np.full(in_range[0], np.nan)
        middle = f(xx_out[in_range])
        end = np.full(xx_out.shape[0] - in_range[-1] - 1, y_in[next_x_in])
        return np.concatenate((front, middle, end))
    elif len(in_range) == 1:
        # Only one point available. Before that point, NaN. Then repeat the
        # single afterwards.
        front = np.full(in_range[0], np.nan)
        end = np.full(xx_out.shape[0] - in_range[0], y_in[next_x_in])
        return np.concatenate((front, end))
    else:
        # No overlap. Either use the most recent value, or NaN.
        if x_in[0] < xx_out[0]:
            return np.full(xx_out.shape[0], y_in[-1])
        else:
            return np.full(xx_out.shape[0], np.nan)


class _LinearInterpolator:
    """An interpolator class. From an earlier iteration of the code and 
    CURRENTLY NOT USED but left for possible future use."""

    def __init__(self):
        self.x_in = None
        self.y_in = None
        self.is_numeric = False

    def calculate(self, x_in, y_in):
        """Do any calculations on the input vectors.

        Args:
            x_in (G3VectorTime) : input x-axis
            y_in (G3Vector or derived) : input y-axis
        """
        if not isinstance(x_in, core.G3VectorTime):
            raise RuntimeError("Argument x_in must be a G3VectorTime instance.")
        self.x_in = x_in
        if not isinstance(y_in, core.G3Vector):
            raise RuntimeError("Argument y_in must be a G3Vector, or " +\
                               "derived class, instance")
        self.y_in = y_in
        for _type in _G3VECTOR_NUMERIC_TYPES:
            if isinstance(self.y_in, _type):
                self.is_numeric = True
                break

    def __call__(self, x_out):
        if x_in == None or y_in == None:
            raise RuntimeError("Interpolation not yet calculated.")

        if not self.is_numeric:
            return _scipy_interpolate(self.x_in, self.y_in, x_out,
                                      kind="previous")


        return _scipy_interpolate(self.x_in, self.y_in, x_out,
                                  kind="linear")

def _full_prov_id(sess, prov_id):
    """Simple book-keeping device: unique handle for session+provider pair."""
    return "%d.%d" % (sess, prov_id)

class _RefTime:
    def __init__(self, t, idx):
        """Simple structure for storing a reference time vector, along with a
        unique index for this reference time chunk."""
        self.t = t
        self.idx = idx

class _Provider:
    def __init__(self, sess, prov_id):
        """A class for tracking block streams of a given provider.

        This class takes in frames from the provider, sends the data to the
        correct block stream, and, when asked, can return frames that have been
        interpolated to the reference timestamps.
        
        Arguments:
            sess : An `HKSessionHelper` object, which we use to create data
              frames.
            prov_id : The provider ID.

        Attributes:
            sess : An `HKSessionHelper`.
            pi : The provider ID
            block_streams : Instances of `_BlockStream` associated with this
              provider.
            have_ref_time : This will get set to True by `has_ref_time()` if the
              reference time vector is in a block stream of this provider.
            ref_time_field : The name of the field in the block stream that the
              user requested provide the reference time vector; set by
              `has_ref_time().
            last_ref_times : If the last call to `add_frame()` had the reference
              time vector, then that time vector is copied to this variable.
              Otherwise it will be set to `None`.
            ref_time_idx_done : A list of reference time indices indicating
              which chunks have been processed.
        """
        self.sess = sess
        self.pi = prov_id
        self.block_streams = {}
        self.have_ref_time = False
        self.ref_time_field = None
        self.last_ref_times = None

    def has_ref_time(self, field):
        """Notify this provider that it has the reference time.
        
        Arguments:
            field : The name of the field connected to the reference time.
        """
        self.have_ref_time = True
        self.ref_time_field = field

    def add_frame(self, f_in):
        """Copy the data from a data frame into the block stream buffers.

        Arguments:
            f_in : The frame.
        """
        # Reset the buffer of new reference times.
        self.last_ref_times = None

        # Process each block in this frame.
        for b, bn in zip(f_in["blocks"], f_in["block_names"]):
            # Add this block to its Block Stream, registering the latter if
            # it does not yet exist.
            full_bn = "%d.%s" % (self.pi, bn)
            if bn not in self.block_streams:
                core.log_info("Creating _BlockStream() %s." % (full_bn))
                self.block_streams[bn] = _BlockStream(bn, self.pi)
            self.block_streams[bn].add_block(b)
                   
            # Check to see if the reference timestamps are in this frame.
            if self.have_ref_time:
                if self.ref_time_field in b.keys():
                    core.log_info("New reference time chunk found.")
                    self.last_ref_times = b.times
                    self.block_streams[bn].have_ref_time = True


    def ref_time_idx_done(self, idx):
        # Check to see if all block streams have been interpolated for this
        # reference time chunk.
        for bs in self.block_streams.values():
            if idx not in bs.ref_time_idx_done:
                return False
        return True

    def process(self, ref_time, flush=False):
        """Compile interpolated frames to be emitted.

        Arguments:
            ref_time : The reference time vector to which all data will be
              interpolated.
            flush : If True, then discard data that cannot be interpolated.

        Returns: A list of frames.
        """
        f_out = None

        # Loop through all the block streams of this provider.
        for bs in self.block_streams.values():
            # Only process the block stream if it has not already been processed
            # for this particular section of the reference time vector.
            if ref_time.idx not in bs.ref_time_idx_done:
                b = bs.process(ref_time, flush=flush)
                if b:
                    # The block stream has emitted some data to be output. Add
                    # it to the frame for emission (creating a frame if it
                    # doesn't yet exist).
                    if not f_out:
                        f_out = self.sess.data_frame(self.pi, ref_time.t[0])
                    f_out["blocks"].append(b)
                    f_out["block_names"].append(bs.name)

        if f_out:
            return [f_out]
        else:
            return []

class _BlockStream:
    def __init__(self, name, prov_id):
        """Buffers data from a block stream that is waiting to be interpolated.

        Arguments:
            name : The name of the block stream (used only for debugging).
            prov_id : The provider ID (used only for debugging).

        Attributes:
            ref_time_idx_done : A list of reference time indices indicating
              which chunks have been processed.
            name : The name of the block stream (used only for debugging).
            pi : The provider ID (used only for debugging).
            have_ref_time : Indicates whether this block stream provides the
              reference time.
        """
        self._data = {}
        self._type = {}
        self.ref_time_idx_done = []
        self.name = name
        self.pi = prov_id
        self._times = np.array([])
        self._data = {}
        self._type = {}
        self.have_ref_time = False

    def add_block(self, block):
        """Add the data from a block into our buffers.

        Arguments:
            block : The block whose data should be added.
        b : A G3TimeSampleMap object.
        """
        self._times = np.append(self._times,
                                np.array([t.time for t in block.times]))
        for k in block.keys():
            if k not in self._data.keys():
                self._data[k] = np.array([])
                self._type[k] = type(block[k])
            self._data[k] = np.append(self._data[k], np.array(block[k]))

    def process(self, ref_time, flush=False):
        """Interpolate data, if possible.

        Arguments:
            ref_time : A `_RefTime` object with the chunk of times over which to
              interpolate.
            flush : If True, emit a block even if the block stream isn't caught
              up yet. (Unless the `ref_time` is empty.)

        Returns: `None` if the block stream isn't caught up to the reference
          time or if no reference time is provided; otherwise, an (interpolated)
          block of data.
        """
        if len(self._times) == 0:
            return None

        # Return nothing if we are not caught up, unless we have been asked to
        # flush.
        ref_t_max = ref_time.t[-1].time
        if self._times[-1] < ref_t_max:
            if not flush:
                 # Should we do the following?
                 #   self.ref_time_idx_done.append(ref_time.idx)
                return None
            else:
                core.log_info("Flushing blockstream %d.%s." %\
                              (self.pi, self.name))

        # Once we are done interpolating, we will be able to discard all data
        # before the index last_i.
        last_i = bisect.bisect_left(self._times, ref_t_max) - 1
        if last_i < 0:
            last_i = 0

        # Make block for output.
        b = core.G3TimesampleMap()
        b.times = core.G3VectorTime(ref_time.t)
            
        if False and self.have_ref_time:
            # TO DO:
            # Special case: this block stream is providing the reference time,
            # so no need for any interpolation.
            i = bisect.bisect_left(self._times, ref_time.t[0].time)
            j = i + len(ref_time.t)
            for k in self._data.keys():
                b[k] = self._type[k](self._data[k][i:j])
                self._data[k] = self._data[k][last_i]
        else:
            for k in self._data.keys():
                # TO DO: test for numeric!
                # Interpolate, stick into the block that will be returned, and
                # then remove the data we no longer need from our buffer.
                b[k] = self._type[k](_scipy_interpolate(self._times,
                                                        self._data[k],
                                                        ref_time.t))
                self._data[k] = self._data[k][last_i:]

        # Remove the timestamps we no longer need from our buffer.
        self._times = self._times[last_i:]

        # Record the fact that we have successfully interpolated this time
        # chunk.
        self.ref_time_idx_done.append(ref_time.idx)

        return b


class HKResampler:
    def __init__(self, t_ref_field):
        """Interpolate HK data frames.

        This is in progress. Currently you can do the following:
        - Interpolate all data so that it has the same time-axis as a field
          chosen by the user (`t_ref_field`)
          + The interpolation is basic: it doesn't try to do anything clever if
            there are large gaps in the data being interpolated.
          + There is not filtering for the case when you are downsampling, so
            there will be antialiasing.
        - It works on the "test_data" from simons1 (barring bug(s) listed
          below).

        Bugs to fix:
        - Currently, if a provider is dropped at the appearance of a new
          status frame, all its data are
          just flushed immediately. However, if there is no recent data from the
          field providing the reference times, then not all the data can be
          processed in the provider that is dropped, and data will be missing.
          The correct solution is to delay writing the status frame (and frames
          from non-dropped providers) until dropped provider can write out all
          its data. But this requires more book-keeping. See comment labelled
          "FIX1".
        - _Possible_ bug: should all providers be flushed when a new session
          starts? Technically frames from the previous session could still
          appear in the stream. See comment labelled "FIX2".

        Functionality to add:
        - In addition to using the times of one of the fields, add the option to
          specify:
          + A vector of times provided by the user.
          + A cadence, with the grid of times automatically generated.
        - Filtering of the data for the case when you are downsampling.
        - More options for how to handle large gaps in the data.

        Efficiencies:
        - Currently everything is done in python with numpy. Can probably get
          speed-ups by writing sections in C++.

        Arguments:
            t_ref_field : The field to use for a reference timeline, in the
              format {full-provider-name}.{field-name}. E.g., 
              observatory.LSA2761.feeds.temperatures.Channel_05_R.
        """
        self._t_ref_prov = ".".join(t_ref_field.split(".")[0:-1])
        self._t_ref_sfield = t_ref_field.split(".")[-1]
        self._sessions = {}
        self._ref_time = []
        self._providers = {}

    def _process_provider(self, prov, flush=False):
        """Process provider ready for interpolation. Return frames for
        output."""
        f_out = []
        for rt in self._ref_time:
            if not prov.ref_time_idx_done(rt.idx):
                f = prov.process(rt, flush=flush)
                # If the provider was not ready to interpolate this time
                # chunk, then do not continue, since we need each block
                # stream to appear in temporal sequence.
                if len(f) == 0 and not flush:
                    break
                f_out += f
        return f_out

    def _process_all_providers(self, flush=False):
        """Processing all providers."""
        core.log_info("Processing providers for new frames.")
        f_out = []
        for prov in self._providers.values():
            f_out += self._process_provider(prov, flush=flush)
        return f_out

    def __call__(self, f_in):
        """Processes a frame. Only Housekeeping frames with key `"hkagg_type"`
        set to `so3g.HKFrameType.data` will be resampled; any other frame will
        be returned as-is.

        Args:
            f: The frame to be processed

        Returns:
            The downsampled frame.
        """
        # If we are at the end: flush everything regardless.
        if f_in.type == core.G3FrameType.EndProcessing:
            f_out = self._process_all_providers(flush=True)
            core.log_info("End of processing: returning %d frames." %\
                          len(f_out))
            return f_out + [f_in]

        # Pass through non-HK frames.
        if f_in.type != core.G3FrameType.Housekeeping:
            return [f_in]

        # Ensure we are using the most recent version.
        hkagg_vers = f_in.get("hkagg_version", 0)
        assert(hkagg_vers == 2)

        # Now figure out what kind of frame we have and act appropriately. In
        # all cases we potentially accumulate frames to be emitted in f_out.
        sess_id = f_in["session_id"]
        f_out = []

        if f_in["hkagg_type"] == so3g.HKFrameType.session:
            # Only update the session and flush our buffers if this is truly
            # a new session and not just the beginning of a new file …
            if sess_id not in self._sessions.keys():
                core.log_info("Registering session frame with ID = %d." %\
                              sess_id)
                hkagg_desc = f_in.get("description", None)
                self._sessions[sess_id] = so3g.hk.HKSessionHelper(sess_id,
                                                     hkagg_version=hkagg_vers,
            
                                                     description=hkagg_desc)
                # Process remaining providers and then pass through the session
                # frame.
                # FIX2? Or not? Perhaps don't flush!
                #f_out += self._process_all_providers(flush=True)
                f_out.append(f_in)
            else:
                core.log_info("Session ID %d already registered. Ignoring "\
                              "frame." % sess_id)

        elif f_in["hkagg_type"] == so3g.HKFrameType.status:
            # Whenever there is a new status frame, go through and register all
            # the providers, and also figure out where our reference timeline 
            # is. Throw an exception if the latter is not present.
            core.log_info("Read status frame: now processing:")
            t_ref_pi = False
            curr_providers = []
            for f in f_in["providers"]:
                # Create our providers object if it doesn't exist.
                pi = f["prov_id"].value
                full_pi = _full_prov_id(sess_id, pi)
                curr_providers.append(pi)
                if full_pi not in self._providers.keys():
                    core.log_info("  Creating _Provider() ID %d (%s)." %\
                                  (pi, f["description"]))
                    self._providers[full_pi] = _Provider(
                                                 self._sessions[sess_id], pi)
                else:
                    core.log_debug("  _Provider() for ID %d (%s) already "\
                                   "exists for this session. Doing nothing." %\
                                   (pi, f["description"]))

                # Check for reference time provider.
                if f["description"].value == self._t_ref_prov:
                    core.log_info("Reference time provider == %s." %\
                                  full_pi)
                    self._providers[full_pi].has_ref_time(self._t_ref_sfield)
                    t_ref_pi = True
            if not t_ref_pi:
                raise RuntimeError("No provider named \"%s\" found." % \
                                   self._t_ref_prov)

            # If any providers have been discontinued, flush them and remove
            # them.
            # FIX1: don't do flush=True. Delay writing out this status frame
            # until we can process the whole provider.
            discontinued = []
            for k, prov in self._providers.items():
                if prov.sess.session_id == sess_id:
                    if prov.pi not in curr_providers:
                        core.log_info("Provider %d discontinued. Flushing." %\
                                      prov.pi)
                        self._process_provider(prov, flush=True)
                        discontinued.append(k)
            for k in discontinued:
                del self._providers[k]

            # Pass the provider frame through.
            # FIX1: delay emitting the provider frame until dropped providers
            # can be finished up.
            f_out.append(f_in)

        elif f_in["hkagg_type"] == so3g.HKFrameType.data:
            # Add the frame data to the relevenant _Provider.
            pi = f_in["prov_id"]
            full_pi = _full_prov_id(sess_id, f_in["prov_id"])
            curr_prov = self._providers[full_pi]
            curr_prov.add_frame(f_in)

            # If we have just gotten a new chunk of reference times, go through
            # and process all our providers to do any new interpolation.
            if curr_prov.have_ref_time:
                if curr_prov.last_ref_times:
                    idx = len(self._ref_time)
                    self._ref_time.append(_RefTime(curr_prov.last_ref_times,
                                                   idx))
                    f_out += self._process_all_providers()
        
        core.log_debug("Returning %d frames." % len(f_out))
        return f_out


if __name__ == '__main__':
    import so3g
    import sys

    # Usage: python resampler.py [output_filename] [input_filenames]…

    core.set_log_level(core.G3LogLevel.LOG_DEBUG)

    path = sys.argv[2:]

    p = core.G3Pipeline()
    p.Add(core.G3Reader(path))
    p.Add(so3g.hk.HKTranslator())
    p.Add(HKResampler("observatory.LSA2761.feeds.temperatures.Channel_05_R"))
#    p.Add(HKResampler("observatory.LSA22YE.feeds.temperatures.Channel_01_R"))
#    p.Add(HKResampler("observatory.bluefors.feeds.bluefors.pressure_ch1"))
    p.Add(core.G3Writer, filename=sys.argv[1])
    p.Run()

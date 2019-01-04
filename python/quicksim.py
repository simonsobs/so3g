"""
Code to help with generating simulated TOD data streams, including
simple signal injection.

Routines are meant to bootstrap off of each other; a schedule block
can be used to induce frames with a scan pattern in them, then the
scan patterns can be passed to the detector data simulator.  Then you
can add realism or whatever.
"""

import so3g
from spt3g import core
from spt3g import coordinateutils as cu3g

import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
#from so3g import coords
from so3g import coords

FT = core.G3FrameType

def create_CES_observation(start_time, stop_time, az_min, az_max, el,
                           scan_speed=1.5):
    """
    Returns a G3Frame(Observation) with a few fields set to describe a
    constant-elevation scan.
    """
    f = core.G3Frame(FT.Observation)
    f['start_time'] = start_time
    f['stop_time'] = stop_time
    f['scan_pattern'] = 'CES'
    f['el'] = el
    f['az_min'] = az_min
    f['az_max'] = az_max
    f['az_vel'] = scan_speed
    return f

def create_focal_plane(n_det):
    # Generate a simple squarish grid.
    i, p = np.arange(n_det) // 2, np.arange(n_det) % 2
    side = max(2, int(i[-1]**.5))
    row, col = i // side, i % side
    pol_fam = (row + col) % 2
    pol = (pol_fam * 45 + p * 90) * core.G3Units.deg
    x = (col/(side-1) - .5) * 1. * core.G3Units.deg
    y = (row/(side-1) - .5) * 1. * core.G3Units.deg

    # Convert to quternions in the prescribed way.
    phi = np.arctan2(y, x)
    theta = np.arcsin((x**2+y**2)**.5 / core.G3Units.rad)
    q = (coords.q_euler(2, phi) *
         coords.q_euler(1, theta) *
         coords.q_euler(2,-phi) *
         coords.q_euler(2, pol / core.G3Units.rad))

    f = core.G3Frame(FT.Calibration)  #booo
    f['cal_type'] = 'focal_plane'
    # For now just store a vector of detector names, then a vector of
    # boresight-relative quaternion rotations for corresponding dets.
    f['signal_q'] = q
    f['signal_x'] = core.G3VectorDouble(x)
    f['signal_y'] = core.G3VectorDouble(y)
    f['signal_theta'] = core.G3VectorDouble(theta)
    f['signal_phi'] = core.G3VectorDouble(phi)
    f['signal_pol'] = core.G3VectorDouble(pol)
    f['signal_names'] = core.G3VectorString()

    for j in range(n_det):
        f['signal_names'].append('det%04i%s' % (i[j], {0: 'A', 1: 'B'}[p[j]]))
    return f

class PipelineSeeder(list):
    """
    A way to introduce statically generated frames into a pipeline.
    Instantiate this as a list of seed Frames, then add it as the
    first Pipeline element.
    """
    def __call__(self, frame_in):
        output = []
        if frame_in is not None:
            output.append(frame_in)
        if len(self):
            output.append(self.pop(0))
        return output

class ScanPatternGenerator:
    """
    Ignores all frames except "Observation".  When it encounters an
    Observation frame, it appends a series of Scan frames containing
    boresight position information.
    """
    def __call__(self, f):
        # Preserve the input frame, even if it's an Observation.
        output = [f]

        if f.type != core.G3FrameType.Observation:
            return output

        assert(f['scan_pattern'] == 'CES')
        freq = 200.  # hz.
        n_samp = int((f['stop_time'] - f['start_time']) * freq)
        time_vec = f['start_time'] + np.arange(n_samp) / freq
        el_vec = np.zeros(n_samp) + f['el']

        swing_time = (f['az_max'] - f['az_min']) / f['az_vel']
        if swing_time < 1.:  # block pathologically quick scans.
            az_vec = np.zeros(n_samp)  + f['az_min']
            swing_points = [0, n_samp]
        else:
            phase_vec = ((time_vec - time_vec[0]) / swing_time + 0.5)
            swing_points = (np.diff(np.floor(phase_vec)) != 0).nonzero()[0]
            swing_points = [0] + list(swing_points) + [n_samp]
            phase_vec = phase_vec % 2
            phase_vec[phase_vec > 1.] = 2. - phase_vec[phase_vec > 1.]
            az_vec = f['az_max'] * (1+phase_vec)/2 + f['az_min'] * (1-phase_vec)/2

        for i0, i1 in zip(swing_points[:-1], swing_points[1:]):
            f0 = core.G3Frame()
            f0.type = core.G3FrameType.Scan
            benc = so3g.IrregBlockDouble()
            benc.t = time_vec[i0:i1]
            benc.data['az'] = az_vec[i0:i1]
            benc.data['el'] = el_vec[i0:i1]
            benc.data['corot'] = el_vec[i0:i1] - 50.
            f0['vertex_enc_raw'] = benc
            output.append(f0)
            
        return output

class SignalInjector:
    """
    Based on a cached focal_plane and some kind of start / end time
    determination, populate the 'signal' map with zeroed timestream
    vectors for all detectors.
    """

    focal_plane = None
    start_time = None

    def __init__(self, f_samp=200.):
        self.f_samp = f_samp
        self.tick_step = max(1, int(np.floor(core.G3Units.sec / f_samp)))
        
    def __call__(self, f):
        if f.type == FT.Calibration and f['cal_type'] == 'focal_plane':
            self.focal_plane = f

        if f.type != FT.Scan:
            return [f]

        # As long as we have a focal_plane, we can create signal vectors.
        if self.focal_plane is None:
            return [f]
        f['signal'] = core.G3TimestreamMap()

        # Determine time samples we will be covering.
        if self.start_time is None:
            first = f['vertex_enc_raw'].t[0] * core.G3Units.sec
            self.start_time = core.G3Time(
                np.ceil(first / self.tick_step) * self.tick_step)
        # And we will end before...
        last = core.G3Time(f['vertex_enc_raw'].t[-1] * core.G3Units.sec)
        n = int((last.time - self.start_time.time) / self.tick_step)
        end_time = core.G3Time(self.start_time.time + n * self.tick_step)

        z = np.zeros(n)
        for k in self.focal_plane['signal_names']:
            f['signal'][k] = core.G3Timestream(z)

        # You can't broadcast-set the start and end time unless the
        # elements are already populated.
        f['signal'].start = self.start_time
        f['signal'].stop = end_time

        self.start_time = end_time
        return [f]

class NaiveBoresightPointingModel:
    """
    Resample raw encoder signal to match detectors (making it a
    primary field).  Define boresight position as a quaternion
    rotation by taking the az and el encoders as true pointing.

    Because the raw encoder signal included with a frame may not
    exactly cover the detector timeline, this module buffers the raw data.
    """
    frame_buffer = None
    raw_buffer = None
    def __init__(self, signal_name='signal',
                 boresight_name='boresight',
                 enc_name='vertex_enc_raw'
    ):
        self.frame_buffer = []
        self.raw_buffer = []
        self.signal_name = signal_name
        self.boresight_name = boresight_name
        self.enc_name = enc_name

    def __call__(self, f):
        if f.type == FT.Calibration and f['cal_type'] == 'focal_plane':
            self.focal_plane = f

        if f.type != FT.Scan:
            return [f]

        if f.type == FT.EndProcessing:
            flush = True
        else:
            flush = False
            self.frame_buffer.append(f)
        self.raw_buffer.append(f[self.enc_name])
        
        # Figure out what frames we're able to process, given info we have.
        frames_out = []

        # Work in units of seconds.
        raw_t0, raw_t1 = self.raw_buffer[0].t[0], self.raw_buffer[-1].t[-1]
        # Process any frame that ends before raw_t1.
        frame_index = 0
        while len(self.frame_buffer) > 0:
            f = self.frame_buffer[0]
            if not flush and (f['signal'].stop.time / core.G3Units.sec > raw_t1):
                break
            sig = f[self.signal_name]  # f['signal']
            frame_t0 = sig.start.time / core.G3Units.sec
            frame_t1 = sig.stop.time / core.G3Units.sec
            tick_rate = sig.sample_rate / core.G3Units.Hz
            # Figure out what range of samples we will be able to set.
            start_index = np.ceil((raw_t0 - frame_t0)*tick_rate)
            end_index = np.floor((raw_t1 - frame_t0)*tick_rate) + 1
            start_index = max(0, int(start_index))
            end_index = min(int(end_index), sig.n_samples)
            if end_index != sig.n_samples and not flush:
                # Buffer.
                break

            # Otherwise, do the interpolations...
            frames_out.append(self.frame_buffer.pop(0))
            t_raw = np.hstack([r.t for r in self.raw_buffer])
            t_int = frame_t0 + np.arange(start_index, end_index) / tick_rate
            boresight = core.G3TimestreamMap()
            vs = {}
            for k in ['az', 'el', 'corot']:
                interp = spline1d(t_raw, np.hstack([r.data[k] for r in self.raw_buffer]))
                v = np.empty(sig.n_samples)
                v[:start_index] = np.nan
                v[start_index:end_index] = interp(t_int)
                v[end_index:] = np.nan
                vs[k] = v
                boresight[k] = core.G3Timestream(vs[k])

            boresight.start = sig.start
            boresight.stop = sig.stop
            f[self.boresight_name] = boresight  # f['boresight']

            # Compute quaternion.
            q = (# Sky <-- near (el, az=0)
                 coords.q_euler(2, -vs['az'] * np.pi/180)  *
                 # ... sky at az=0 <-- near (el=0,az=0)
                 coords.q_euler(1, -vs['el'] * np.pi/180) *
                 # ... (1,-xi,eta) <-- (-eta,-xi,1)
                 coords.q_euler(1, np.pi/2) *
                 # ... (-eta,-xi,1) <-- (eta,xi,1)
                 coords.q_euler(2, np.pi)
            )
            # Note that there's no "TimestreamQuat" class.  So no timestamps.
            f[self.boresight_name + '_q'] = q   # f['boresight_q']
            
        # Discard raw data we're not using any more.  Out of caution,
        # keep one more frame than we have buffered.
        while len(self.raw_buffer)  - len(self.frame_buffer) > 2:
            self.raw_buffer.pop(0)

        return frames_out

class Inspector:
    def __call__(self, f):
        print('Inspector: %s' % f)
        print('The frame is called "f".')
        import pdb
        pdb.set_trace()
        
if __name__ == '__main__':
    #core.set_log_level(core.G3LogLevel.LOG_TRACE)

    test_file = 'sim_out.g3'
    print('Streaming to %s' % test_file)

    start_time = core.G3Time('2022-06-01T00:00:00')  #isoformat
    start_ctime = start_time.time / core.G3Units.sec
    length_s = 60 * 15

    p = core.G3Pipeline()

    p.Add(PipelineSeeder([
        create_focal_plane(200),
        create_CES_observation(start_ctime,start_ctime+length_s,
                               45,66,45),
    ]))
    
    p.Add(ScanPatternGenerator)
    p.Add(SignalInjector)
    #p.Add(Inspector)
    p.Add(NaiveBoresightPointingModel)

    p.Add(core.G3Writer, filename=test_file)

    p.Run()
    del p
    
    print('Reading back:')
    for f in core.G3File(test_file):
        if f.type == FT.Observation:
            print(f, f['start_time'])
        if f.type == FT.Scan:
            print(f, f['vertex_enc_raw'].t[0])

    print('Done reading %s.' % test_file)

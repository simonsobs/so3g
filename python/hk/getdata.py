"""Interface for querying and loading contiguous vectors from
Housekeeping frame streams.  We want to expose a sisock-compatible
API, where time ranges are an essential part of every query.

Use HKArchiveScanner to scan a set of G3 files containing Housekeeping
frames.  Use the resulting HKArchive to perform random I/O over a
sisock-style API.

"""
import itertools
import os
import pytz
import yaml
import logging
import pickle
import glob

import numpy as np
import datetime as dt


import so3g
from spt3g import core


hk_logger = logging.getLogger(__name__)
hk_logger.setLevel(logging.INFO)

SPAN_BUFFER_SECONDS = 1.0

def is_sub_seq(full_seq, sub_seq):
    """Return true if sub_seq is a sub-sequence of full_seq.

    """
    if len(sub_seq) == 0:
        return True
    for idx0,x0 in enumerate(full_seq):
        if x0 == sub_seq[0]:
            if is_sub_seq(full_seq[idx0+1:], sub_seq[1:]):
                return True
    return False


_SCHEMA_V1_BLOCK_TYPES = [
    core.G3VectorDouble,
    core.G3VectorInt,
    core.G3VectorString,
    core.G3VectorBool,
]


def _concat_hk_stream(blocks_in):
    """Concatenates an ordered list of compatible HK blocks into a single
    frame.  Each block should be a valid G3TimesampleMap with the same
    keys.

    Returns a single G3TimesampleMap with all fields concatenated.

    """
    blk = core.G3TimesampleMap()
    blk.times = core.G3VectorTime(blocks_in[0].times)
    fields = list(blocks_in[0].keys())
    for f in fields:
        f_short = f.split('.')[-1]
        blk[f] = blocks_in[0][f_short]
    for b in blocks_in[1:]:
        blk.times.extend(b.times)
    for f in fields:
        f_short = f.split('.')[-1]
        for _type in _SCHEMA_V1_BLOCK_TYPES:
            if isinstance(blocks_in[0][f_short], _type):
                break
        else:
            raise RuntimeError('Field "%s" is of unsupported type %s.' %
                               (f_short, type(blocks_in[0][f_short])))
        for b in blocks_in[1:]:
            blk[f].extend(b[f_short])
    return blk


class HKArchive:
    """Container for information necessary to determine what data fields
    are present in a data archive at what times.  This object has
    methods that can determine what fields have data over a given time
    range, and can group fields that share a timeline (i.e. are
    co-sampled) over that range.

    """
    field_groups = None
    def __init__(self, field_groups=None):
        if field_groups is not None:
            self.field_groups = list(field_groups)
        # A translator is used to update frames, on the fly, to the
        # modern schema assumed here.
        self.translator = so3g.hk.HKTranslator()

    def _get_groups(self, fields=None, start=None, end=None,
                    short_match=False):
        """Helper for get_fields and get_data.  Determines which fields, of
        those listed in fields, are present in the archive between the
        specified times.

        Args:
            fields (list of strings): List of field names.  If None,
              all fields are considered.
            start (timestamp): Earliest time to consider.  If None,
              consider arbitrarily early times.
            end (timestamp): Latest time to consider.  If None,
              consider arbitrarily late times.  Note that ``start``
              and ``end`` form a semi-closed interval that includes
              start but excludes end.
            short_match (bool): If True, then a requested field will
              be considered to match an archive field if its
              "."-tokenized form is a sub-sequence of the
              "."-tokenized field in the archive.  For example, the
              archive field "obs.agentX.feeds.therms.c1" may be
              matched by the requested field "agentX.c1".  In the case
              that multiple archive fields match a requested field, a
              ValueError is thrown.

        Returns:
            List of (group_name, group_fields, fgs).  The
            ``group_name`` is internally generated... probably
            something like "group0".  ``group_fields`` is a list of
            fields belonging to this group.  ``fgs`` is a list of
            _FieldGroup objects relevant to the fields in this
            group.

        Notes:
            Each entry in the returned list carries information for
            set of fields that are co-sampled over the requested time
            interval.  Each of the requested fields will thus appear
            in at most one group.

            If a requested field is not found, it will not be listed
            in any group.

        """
        span = so3g.IntervalsDouble()
        if start is None:
            start = span.domain[0]
        if end is None:
            end = span.domain[1]
        span.add_interval(start, end)

        if short_match and (fields is not None):
            field_seqs = [f.split('.') for f in fields]
            short_match_map = {}  # map from shortened name to full field name.

        field_map = {}
        for fg in self.field_groups:
            both = span * fg.cover
            if len(both.array()) > 0:
                for f in fg.fields:
                    full_field = fg.prefix + '.' + f
                    key_field = full_field
                    if fields is not None:
                        # User is interested only in particular fields.
                        if short_match:
                            for seq in field_seqs:
                                if is_sub_seq(full_field.split('.'), seq):
                                    key_field = '.'.join(seq)
                                    prev_short_match = short_match_map.get(key_field)
                                    if prev_short_match not in [None, full_field]:
                                        raise ValueError("Multiple matches for soft key: %s [%s, %s]." %
                                                         (key_field, prev_short_match, full_field))
                                    short_match_map[key_field] = full_field
                                    break
                            else:
                                continue
                        elif full_field not in fields:
                            continue
                    if not key_field in field_map:
                        field_map[key_field] = [fg]
                    else:
                        field_map[key_field].append(fg)

        # Sort each list of field_groups by object id -- all we care
        # about is whether two fields have the same field group set.
        [f.sort(key=lambda obj: id(obj)) for f in field_map.values()]

        # Now group together fields if they have identical
        # field_group lists (because when they do, they can share a
        # timeline).
        keys = list(field_map.keys())
        grouped = []
        i0, i1 = 0, 1
        while i0 < len(keys):
            while i1 < len(keys) and field_map[keys[i0]] == field_map[keys[i1]]:
                i1 += 1
            group_keys = keys[i0:i1]
            group_name = 'group%i' % len(grouped)
            grouped.append((group_name, group_keys, field_map[group_keys[0]]))
            i0, i1 = i1, i1+1
        return grouped

    def get_fields(self, start=None, end=None):
        """Get list of fields that might have a sample in the time interval
        [start,end).

        Returns the pair of dictionaries ``(fields, timelines)``.

        The ``fields`` dictionary is a map from field name to a block
        of field information.  The ``timelines`` dictionary is a map
        from timeline name to a block of timeline information.

        """
        grouped = self._get_groups(None, start, end)
        # Construct the fields and timelines dicts.
        field_map, timeline_map = {}, {}
        for (group_name, fields, data) in grouped:
            timeline_map[group_name] = {'interval': None, 'field': fields}
            for f in fields:
                field_map[f]  = {
                    'description': None,
                    'timeline': group_name,
                    'type': 'number',
                    'units': None,
                }
        return field_map, timeline_map

    def get_data(self, field=None, start=None, end=None, min_stride=None,
                 raw=False, short_match=False):
        """Load data from specified field(s) between specified times.

        Arguments ``field``, ``start``, ``end``, ``short_match`` are
        as described in _get_groups.

        Arguments:
            min_stride (float): Specifies the minimum sample spacing,
                in seconds.  Ignored in this implementation.
            raw (bool): If true, return G3 blocks instead of numpy
                arrays.

        Returns:
            Pair of dictionaries, (data, timelines).  The ``data``
            dictionary is a simple map from field name to a numpy
            array of readings.  The ``timelines`` dictionary is a map
            from field group name to a dictionary of timeline
            information, which has entries:

            - ``'t'``: numpy array of timestamps
            - ``'fields'``: list of fields belonging to this group.
            - ``'finalized_until'``: in cases where the data are still
              in flux, this field provides a timestamp that may be
              taken as the earliest time that needs to be requeried.
              This is part of the interface in order to support data
              streams that are being updated in real time.

            If user requested raw=True, then return value is a list
            of tuples of the form (group_name, block) where block is a
            single G3TimesampleMap carrying all the data for that
            co-sampled group.

        """
        grouped = self._get_groups(field, start, end, short_match=short_match)
        hk_logger.debug('get_data: _get_groups returns %i groups.' % len(grouped))

        # Pass through the metadata.  Collect information on what
        # field groups are present in what frames of what files; sort
        # that info by file and offset so we make a single monotonic
        # pass through the frames.
        group_info = {
            #  group_name: {'types': [dtype, ...],
            #               'fields': [(full_name, short_name), ...],
            #               'count': n},
            #  ...
        }
        files = {
            # filename: {
            #   offset: [(block_index, group_name, output_offset), ...]],
            #   ...
            #   },
            # ...
        }

        def check_overlap(time_range):
            # Note the time_range is inclusive on both ends.
            return ((start is None or start <= time_range[1]) and
                    (end is None or end > time_range[0]))

        for group_name, fields, fgrps in grouped:
            # This is a group of co-sampled fields.  The fields share
            # a sample count and a frame-index map.
            all_frame_refs = []
            for fg in fgrps:
                all_frame_refs.extend(
                    [(b['time_range'], b['count'], b['filename'], b['byte_offset'], b['block_index'])
                     for b in fg.index_info])
            all_frame_refs.sort()
            vector_offset = 0
            for time_range, n, filename, byte_offset, block_index in all_frame_refs:
                if not check_overlap(time_range):
                    continue
                if not filename in files:
                    files[filename] = {}
                if byte_offset not in files[filename]:
                    files[filename][byte_offset] = []
                files[filename][byte_offset].append((block_index, group_name, vector_offset))
                vector_offset += n
            group_info[group_name] = {
                'count': vector_offset,
                'fields': [(f, f.split('.')[-1]) for f in fields],
                'types': [],
            }

        # Pass through the data.  Should read the files in order,
        # jumping monotonically through the needed frames.  The data
        # type of output arrays is determined from whatever
        # np.array(G3Object) returns on the first block.  Note strings
        # ('U') have to be handled differently because we can't know
        # the maximum string length from the first block.
        data = {}
        timelines = {}
        for filename, file_map in sorted(files.items()):
            hk_logger.debug('get_data: reading %s' % filename)
            reader = core.G3Reader(filename)
            for byte_offset, frame_info in sorted(file_map.items()):
                # Seek and decode.
                hk_logger.debug('get_data: seeking to %i for %i block extractions' %
                                (byte_offset, len(frame_info)))
                reader.seek(byte_offset)
                frames = reader.Process(None)
                assert(len(frames) == 1)
                frames = self.translator(frames[0])
                frame = frames[0]
                # Now expand all blocks.
                for block_index, group_name, offset in frame_info:
                    block = frame['blocks'][block_index]
                    gi = group_info[group_name]
                    if raw:
                        # Short-circuit if in raw mode, just collect all blocks.
                        if group_name not in data:
                            data[group_name] = {}
                            for field, f_short in gi['fields']:
                                data[group_name] = []
                        data[group_name].append(block)
                        continue
                    if group_name not in timelines:
                        # This block is init that happens only once per group.
                        timelines[group_name] = {
                            't_g3': np.zeros(gi['count'], dtype=np.int64),
                            'fields': [f for f,s in gi['fields']],
                        }
                        hk_logger.debug('get_data: creating group "%s" with %i fields'
                                        % (group_name, len(gi['fields'])))
                        # Determine data type of each field and create output arrays.
                        for field, f_short in gi['fields']:
                            dtype = np.array(block[f_short]).dtype
                            gi['types'].append(dtype)
                            if dtype.char == 'U':
                                data[field] = []
                            else:
                                data[field] = np.empty(gi['count'], dtype)
                            hk_logger.debug('get_data: field "%s" has type %s' % (
                                f_short, dtype))
                    # Copy in block data.
                    n = len(block.times)
                    # Note this is in G3 time units for now... fixed later.
                    timelines[group_name]['t_g3'][offset:offset+n] = block.times
                    for (field, f_short), dtype in zip(gi['fields'], gi['types']):
                        if dtype.char == 'U':
                            data[field].append((offset, list(map(str, block[f_short]))))
                        else:
                            # This is a relatively quick copy because
                            # of buffer pass-through from G3... don't
                            # hit the RHS with np.array!
                            data[field][offset:offset+n] = block[f_short]

        if raw:
            return [(group_name, _concat_hk_stream(data[group_name]))
                    for group_name, _, _ in grouped]

        # Finalize string fields.
        for group_name, fields, fgrps in grouped:
            gi = group_info[group_name]
            for (field, f_short), dtype in zip(gi['fields'], gi['types']):
                if dtype.char == 'U':
                    data[field] = np.array(list(itertools.chain(
                        *[x for i, x in sorted(data[field])])))
                    assert(len(data[field]) == gi['count'])

        # Scale out time units.
        for timeline in timelines.values():
            timeline['t'] = timeline.pop('t_g3') / core.G3Units.seconds

        # Restrict to only the requested time range.
        if start is not None or end is not None:
            for timeline in timelines.values():
                i0, i1 = 0, len(timeline['t'])
                if start is not None:
                    i0 = np.searchsorted(timeline['t'], start)
                if end is not None:
                    i1 = np.searchsorted(timeline['t'], end)
                sl = slice(i0, i1)
                timeline['t'] = timeline['t'][sl]
                for k in timeline['fields']:
                    data[k] = data[k][sl]

        # Mark last time
        for timeline in timelines.values():
            if len(timeline['t']):
                timeline['finalized_until'] = timeline['t'][-1]
            else:
                timeline['finalized_until'] = start if start is not None else 0.

        return (data, timelines)

    def simple(self, fields=None, start=None, end=None, min_stride=None,
               raw=False, short_match=True):
        """Load data from specified field(s) between specified times, and
        unpack the data for ease of use.  Use get_data if you want to
        preserve the co-sampling structure.

        Arguments ``field``, ``start``, ``end``, ``short_match`` are
        as described in _get_groups.  However, ``fields`` can be a
        single string rather than a list of strings.

        Note that ``short_match`` defaults to True (which is not the
        case for getdata).x

        Returns:
            List of pairs of numpy arrays (t,y) corresponding to each
            field in the ``fields`` list.  If ``fields`` is a string,
            a simple pair (t,y) is returned.  ``t`` and ``y`` are
            numpy arrays of equal length containing the timestamps and
            field readings, respectively.  In cases where two fields
            are co-sampled, the time vector will be the same object.

            In cases where there are no data for the requested field
            in the time range, a pair of length 0 float arrays is returned.

        """
        unpack = isinstance(fields, str)
        if unpack:
            fields = [fields]
        data, timelines = self.get_data(fields, start, end, min_stride, raw, short_match)
        output = {}
        fields_not_found = [f for f in fields]
        for timeline in timelines.values():
            # Make the array here, so that the same object is returned
            # for all fields in this group.
            _t = np.array(timeline['t'])
            for f in timeline['fields']:
                output[f] = (_t, np.array(data[f]))
                fields_not_found.remove(f)
        nothing = np.zeros((0,))
        for f in fields_not_found:
            output[f] = (nothing, nothing)
        output = [output[f] for f in fields]
        if unpack:
            output = output[0]
        return output


class _HKProvider:
    def __init__(self, prov_id, prefix):
        self.prov_id = prov_id
        self.prefix = prefix
        self.blocks = {}

    @classmethod
    def from_g3(cls, element):
        prov_id = element['prov_id'].value
        prefix = element['description'].value
        return cls(prov_id, prefix)


class HKArchiveScanner:
    """Consumes SO Housekeeping streams and creates a record of what fields
    cover what time ranges.  This can run as a G3Pipeline module, but
    will not be able to record stream indexing information in that
    case.  If it's populated through the process_file method, then
    index information (in the sense of filenames and byte offsets)
    will be stored.

    After processing frames, calling .finalize() will return an
    HKArchive that can be used to load data more efficiently.

    """
    def __init__(self, pre_proc_dir=None, pre_proc_mode=None):
        self.session_id = None
        self.providers = {}
        self.field_groups = []
        self.frame_info = []
        self.counter = -1
        self.translator = so3g.hk.HKTranslator()
        self.pre_proc_dir = pre_proc_dir
        self.pre_proc_mode = pre_proc_mode

    def __call__(self, *args, **kw):
        return self.Process(*args, **kw)

    def Process(self, f, index_info=None):
        """Processes a frame.  Only Housekeeping frames will be examined;
        other frames will simply be counted.  All frames are passed
        through unmodified.  The index_info will be stored along with
        a description of the frame's data; see the .process_file
        function.

        """
        self.counter += 1
        if index_info is None:
            index_info = {'counter': self.counter}

        f = self.translator(f)
        assert(len(f) == 1)
        f = f[0]

        if f.type == core.G3FrameType.EndProcessing:
            return [f]

        if f.type != core.G3FrameType.Housekeeping:
            return [f]

        vers = f.get('hkagg_version', 0)
        assert(vers == 2)

        if f['hkagg_type'] == so3g.HKFrameType.session:
            session_id = f['session_id']
            if self.session_id is not None:
                if self.session_id != session_id:
                    self.flush()  # Note that this sets self.session_id=None.
            if self.session_id is None:
                core.log_info('New HK Session id = %i, timestamp = %i' %
                              (session_id, f['start_time']), unit='HKScanner')
                self.session_id = session_id

        elif f['hkagg_type'] == so3g.HKFrameType.status:
            # If a provider has disappeared, flush its information into a
            # FieldGroup.
            prov_cands = [_HKProvider.from_g3(p) for p in f['providers']]
            to_flush = list(self.providers.keys())  # prov_ids...
            for p in prov_cands:
                if p.prov_id in to_flush:
                    to_flush.remove(p.prov_id) # no, don't.
                else:
                    self.providers[p.prov_id] = p
            for prov_id in to_flush:
                self.flush([prov_id])

        elif f['hkagg_type'] == so3g.HKFrameType.data:
            # Data frame -- merge info for this provider.
            prov = self.providers[f['prov_id']]
            representatives = prov.blocks.keys()

            for bidx, (bname, b) in enumerate(zip(f['block_names'], f['blocks'])):
                assert(isinstance(b, core.G3TimesampleMap))
                if bname not in prov.blocks:
                    prov.blocks[bname] = {'fields': list(b.keys()),
                                          'start': b.times[0].time / core.G3Units.seconds,
                                          'index_info': []}
                # To ensure that the last sample is actually included
                # in the semi-open intervals we use to track frames,
                # the "end" time has to be after the final sample.
                prov.blocks[bname]['end'] = b.times[-1].time / core.G3Units.seconds + SPAN_BUFFER_SECONDS
                ii = {'block_index': bidx,
                      'time_range': (b.times[0].time / core.G3Units.seconds,
                                     b.times[-1].time / core.G3Units.seconds),
                      'count': len(b.times)}
                ii.update(index_info)
                prov.blocks[bname]['index_info'].append(ii)
                
        else:
            core.log_warn('Weird hkagg_type: %i' % f['hkagg_type'],
                          unit='HKScanner')
        return [f]

    def flush(self, provs=None):
        """Convert any cached provider information into _FieldGroup
        information.  Delete the provider information.  This will be
        called automatically during frame processing if a provider
        session ends.  Once frame processing has completed, this
        function should be called, with no arguments, to flush any
        remaining provider sessions.

        Args:
            provs (list or None): list of provider identifiers to
              flush.  If None, flush all and also reset the
              self.session_id.

        Returns:
            None

        """
        if provs is None:
            provs = list(self.providers.keys())
            self.session_id = None
        for p in provs:
            prov = self.providers.pop(p)
            blocks = prov.blocks
            for k, info in blocks.items():
                fg = _FieldGroup(prov.prefix, info['fields'], info['start'],
                                 info['end'], info['index_info'])
                self.field_groups.append(fg)

    def finalize(self):
        """Finalize the processing by calling flush(), then return an
        HKArchive with all the _FieldGroup information from this
        scanning session.
        """
        self.flush()
        return HKArchive(self.field_groups)

    def process_file(self, filename, flush_after=True):
        """Process the file specified by ``filename`` using a G3Reader.
        Each frame from the file is passed to self.Process, with the
        optional index_info argument set to a dictionary containing
        the filename and byte_offset of the frame.

        Internal data grouping will be somewhat cleaner if the
        multiple files from a single aggregator "session" are passed
        to this function in acquisition order.  In that case, call
        with flush_after=False.

        """
        reader = core.G3Reader(filename)
        while True:
            info = {'filename': filename,
                    'byte_offset': reader.tell()}
            frames = reader.Process(None)
            assert len(frames) <= 1
            if len(frames) == 0:
                break
            self.Process(frames[0], info)
        # Calling flush() here protects us against the odd case that
        # we process files from a single session in non-consecutive
        # order.  In that case the start' and 'end' times will get
        # messed up because we can't tell the stream has been
        # re-initialized.
        if flush_after:
            self.flush()


    def process_file_with_cache(self, filename):
        """Processes file specified by ``filename`` using the process_file
           method above. If self.pre_proc_dir is specified (not None), it
           will load pickled HKArchiveScanner objects and concatenates with
           self instead of re-processing each frame, if the corresponding
           file exists.  If the pkl file does not exist, it processes it and
           saves the result (in the pre_proc_dir) so it can be used in the
           future.  If self.pre_proc_dir is not specified, this becomes
           equivalent to process_file.
        """
        if self.pre_proc_dir is None:
            self.process_file(filename)
            return

        folder = os.path.basename(filename)[:5]
        path = os.path.join( self.pre_proc_dir, folder, os.path.basename(filename).replace(".g3",'.pkl') )

        if os.path.exists(path):
            with open(path, 'rb') as pkfl:
                hksc = pickle.load(pkfl)

        else:
            hksc = HKArchiveScanner()
            hksc.process_file(filename)
            # Make dirs if needed
            if not os.path.exists( os.path.dirname(path) ):
                os.makedirs( os.path.dirname(path) )
                if self.pre_proc_mode is not None:
                    os.chmod( os.path.dirname(path), self.pre_proc_mode )
            # Save pkl file
            with open(path, 'wb') as pkfl:
                pickle.dump(hksc, pkfl)
            if self.pre_proc_mode is not None:
                os.chmod( path, self.pre_proc_mode )            

        self.field_groups += hksc.field_groups
        self.counter += hksc.counter




class _FieldGroup:
    """Container object for look-up information associated with a group of
    fields that share a timeline (i.e. a group of fields that are
    co-sampled over some requested time range).

    Attributes:
        prefix (str): Not used.
        fields (list of str): The field names.
        cover (IntervalsDouble): The time range (as unix timestamps)
          over which these fields have data.  This is expressed as an
          IntervalsDouble to enable the use of Intervals arithmetic.
          The range is created from the "start" and "end" arguments of
          the constructor.
        index_info (list of dict): Information that the consumer will
          use to locate and load the data efficiently.  The entries in
          the list represent time-ordered frames. See Notes.

    Notes:

      Each entry of index_info is a dict providing information about a
      single frame and block where data can be found.  The fields are:

      - 'filename' (str): The file in which the frame is located.
      - 'byte_offset' (int): The offset within the file at which to
        find the frame.
      - 'block_index' (int): The index of the block, within the
        corresponding frame, where the data are found.
      - 'count' (int): the number of samples in this block.
      - 'time_range' (tuple): (first time, last time), as unix
        timestamps.  Note the last time is the time of the last
        sample, not some time beyond that.
      - 'counter' (int): An index providing the order in which the
        frames were scanned.

    """
    def __init__(self, prefix, fields, start, end, index_info):
        self.prefix = prefix
        self.fields = list(fields)
        self.cover = so3g.IntervalsDouble().add_interval(start, end)
        self.index_info = index_info
    def __repr__(self):
        try:
            return '_FieldGroup("%s", %i fields, %i segments)' % (
                self.prefix, len(self.fields), len(self.index_info))
        except:
            return '_FieldGroup(<bad internal state!>)'


def to_timestamp(some_time, str_format=None): 
    """Convert the argument to a unix timestamp.

    Args:
      some_time: If a datetime, it is converted to UTC timezone and
        then to a unix timestamp.  If int or float, the value is
        returned unprocessed.  If str, a date will be extracted based
        on a few trial date format strings.
      str_format: a format string (for strptime) to try, instead of
        the default(s).

    Returns:
        float: Unix timestamp corresponding to some_time.

    """
    
    if type(some_time) == dt.datetime:
        return some_time.astimezone(dt.timezone.utc).timestamp()
    if type(some_time) == int or type(some_time) == float:
        return some_time
    if type(some_time) == str:
        if str_format is not None:
            return to_timestamp(pytz.utc.localize( dt.datetime.strptime(some_time, str_format )))
        str_options = ['%Y-%m-%d', '%Y-%m-%d %H:%M', '%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M:%S.%f']
        for option in str_options:
            try:
                return to_timestamp(pytz.utc.localize( dt.datetime.strptime(some_time, option )))
            except:
                continue
        raise ValueError('Could not process string into date object, options are: {}'.format(str_options))
        
    raise ValueError('Type of date / time indication is invalid, accepts datetime, int, float, and string')

def load_range(start, stop, fields=None, alias=None, 
               data_dir=None, config=None, pre_proc_dir=None, pre_proc_mode=None,
               folder_patterns=None, strict=True):
    """Args:

      start: Earliest time to search for data (see note on time
        formats).
      stop: Latest time to search for data (see note on time formats).
      fields: Fields to return, if None, returns all fields.
      alias: If not None, must be a list of strings providing exactly
        one value for each entry in fields.
      data_dir: directory where all the ctime folders are.  If None,
        tries to use $OCS_DATA_DIR.
      config: filename of a .yaml file for loading data_dir / fields /
        alias
      pre_proc_dir: Place to store pickled HKArchiveScanners for g3
        files to speed up loading
      pre_proc_mode: Permissions (passed to os.chmod) to be used on
        dirs and pkl files in the pre_proc_dir. No chmod if None.
      folder_patterns:  List of patterns to search for in folders. If
        None, default pattern is ['{folder}', 'hk_{folder}_*']. If not
        None, usage for .g3 folder: ['{folder}'], and example usage
        for HK books: ['hk_{folder}_lat'] where {folder} will be replaced
        with the first 5 digits of the unix timestamp when looping through
        files.
      strict: If False, log and skip missing fields rather than
        raising a KeyError.

    Returns:

      Dictionary with structure::

        {
            alias[i] : (time[i], data[i])
        }

      It will be masked to only have data between start and stop.

    Notes:

      The "start" and "stop" argument accept a variety of formats,
      including datetime objects, unix timestamps, and strings (see
      to_timestamp function).  In the case of datetime objects, you
      should set tzinfo=dt.timezone.utc explicitly if the system is
      not set to UTC time.

      Example usage::

        fields = [
            'observatory.HAN.feeds.temperatures.Channel 1 T',
            'observatory.HAN.feeds.temperatures.Channel 2 T',
        ]

        alias = [
            'HAN 1', 'HAN 2',
        ]

        folder_patterns = ['hk_{folder}_satp3']

        start = dt.datetime(2020,2,19,18,48)
        stop = dt.datetime(2020,2,22)
        data = load_range(start, stop, fields=fields, alias=alias, folder_patterns=folder_patterns)

        plt.figure()
        for name in alias:
            plt.plot( data[name][0], data[name][1])

    """
    if config is not None:
        if not (data_dir is None and fields is None and alias is None):
            hk_logger.warning('''load_range has a config file - data_dir, fields, and alias are ignored''')
        with open(config, 'r') as f:
            setup = yaml.load(f, Loader=yaml.FullLoader)
        
        if 'data_dir' not in setup.keys():
            raise ValueError('load_range config file requires data_dir entry')
        data_dir = setup['data_dir']
        if 'field_list' not in setup.keys():
            raise ValueError('load_range config file requires field_list entry')
        fields = []
        alias = []
        for k in setup['field_list']:
            fields.append( setup['field_list'][k])
            alias.append( k )
            
    if data_dir is None and 'OCS_DATA_DIR' not in os.environ.keys():
        raise ValueError('if $OCS_DATA_DIR is not defined a data directory must be passed to getdata')
    if data_dir is None:
        data_dir = os.environ['OCS_DATA_DIR']

    hk_logger.debug('Loading data from {}'.format(data_dir))
    
    start_ctime = to_timestamp(start) - 3600
    stop_ctime = to_timestamp(stop) + 3600

    hksc = HKArchiveScanner(pre_proc_dir=pre_proc_dir)


    if folder_patterns is None:
        folder_patterns = ['{folder}', 'hk_{folder}_*']
    for folder in range( int(start_ctime/1e5), int(stop_ctime/1e5)+1):
        bases = []
        for pattern in folder_patterns:
            extended_pattern = pattern.format(folder=folder)

            base = glob.glob(os.path.join(data_dir, extended_pattern))
            bases.extend(base)

        if len(bases) > 1:
            bases.sort
            base = bases[0]
            hk_logger.warn(f"Multiple base folders were found for {folder}. The first one, alphabetically, is selected: {base}")
        elif len(bases) == 1:
            base = bases[0]
        elif len(bases) == 0:
            hk_logger.debug(f"No base folder found for {folder}, skipping")
            continue

        for file in sorted(os.listdir(base)):
            if file.endswith('.yaml'):
                continue
            try:
                t = int(file[:-3])
            except:
                hk_logger.debug('{} does not have the right format, skipping'.format(file))
                continue
            if t >= start_ctime-3600 and t <=stop_ctime+3600:
                hk_logger.debug('Processing {}'.format(base+'/'+file))
                hksc.process_file_with_cache( base+'/'+file)

    
    cat = hksc.finalize()
    start_ctime = to_timestamp(start)
    stop_ctime = to_timestamp(stop)
    
    all_fields,_ = cat.get_fields()
    
    if fields is None:
        fields = all_fields
    if alias is not None:
        if len(alias) != len(fields):
            hk_logger.error('if provided, alias needs to be the length of fields')
    else:
        alias = fields
    
    # Single pass load.
    keepers = []
    for name, field in zip(alias, fields):
        if field not in all_fields:
            hk_logger.debug('`{}` not in available fields, skipping'.format(field))
            continue
        keepers.append((name, field))
    data = cat.simple([f for n, f in keepers],
                      start=start_ctime, end=stop_ctime)
    data = {name: data[i] for i, (name, field) in enumerate(keepers)}

    return data

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        usage='This demo can be used to plot data from SO HK files.')
    parser.add_argument('--field')
    parser.add_argument('files', nargs='+', help=
                        "SO Housekeeping files to load.")
    args = parser.parse_args()

    # Scan the file set.
    hkas = HKArchiveScanner()
    for filename in args.files:
        print(f'Scanning {filename}...')
        hkas.process_file(filename)

    print(' ... finalizing.')
    cat = hkas.finalize()

    # Get list of fields, timelines, spanning all times:
    fields, timelines = cat.get_fields()

    # Show
    for field_name in sorted(fields):
        group_name = fields[field_name]['timeline']
        print(field_name, timelines[group_name]['interval'])

    if args.field is None:
        full_name = list(fields.keys())[0]
        print('Pretending interest in:', full_name)
    else:
        full_name = args.field
        print('User requests:', full_name)

    # Identify the shortest form of this field that also works.
    f_toks = full_name.split('.')
    field_name = full_name
    for i in range(1, 2**len(f_toks)):
        short_name = '.'.join([t for b,t in enumerate(f_toks) if (i >> b) & 1])
        try:
            fields, timelines = cat.get_data([short_name], short_match=True)
        except Exception:
            continue
        if len(short_name) < len(field_name):
            field_name = short_name
            print(field_name)

    print('Name shortened to:', field_name)

    # This is the awkward way, which preserves co-sampled
    # relationships (and is thus annoying to decode in simple cases).
    fields, timelines = cat.get_data([field_name], short_match=True)
    x0 = list(timelines.values())[0]['t']
    y0 = fields[field_name]

    # This is the easy way, which just gives you one timeline per
    # requested field.
    x1, y1 = cat.simple(field_name)
    
    assert np.all(np.array(x0) == x1) and np.all(np.array(y0) == y1)

    import pylab as pl
    pl.plot(x1, y1)
    pl.title(field_name)
    pl.show()

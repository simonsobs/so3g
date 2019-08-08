"""Interface for querying and loading contiguous vectors from
Housekeeping frame streams.  We want to expose a sisock-compatible
API, where time ranges are an essential part of every query.

Use HKArchiveScanner to scan a set of G3 files containing Housekeeping
frames.  Use the resulting HKArchive to perform random I/O over a
sisock-style API.

"""

import so3g
from spt3g import core
import numpy as np

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


class HKArchive:

    """Contains information necessary to determine what data fields are
    present in a data archive at what times.  It also knows how to
    group fields that share a commong timeline.
    """
    field_groups = None
    def __init__(self, field_groups=None):
        if field_groups is not None:
            self.field_groups = list(field_groups)

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

        """
        grouped = self._get_groups(field, start, end, short_match=short_match)
        handles = {}  # filename -> G3IndexedReader map.
        blocks_out = []
        for group_name, fields, fgrps in grouped:
            blocks_in = []
            for fg in fgrps:
                for r in fg.index_info:
                    fn,off = r['filename'], r['byte_offset']
                    if not fn in handles:
                        handles[fn] = so3g.G3IndexedReader(fn)
                    handles[fn].Seek(off)
                    fn = handles[fn].Process(None)
                    assert(len(fn) == 1)
                    # Find the right block.
                    for blk in fn[0]['blocks']:
                        test_f = fields[0].split('.')[-1]   ## dump prefix.
                        if test_f in blk.data.keys():
                            blocks_in.append(blk)
                            break
            # Sort those blocks by timestamp. (Otherwise they'll stay sorted by object id :)
            blocks_in.sort(key=lambda b: b.t[0])
            # Create a new Block for this group.
            blk = so3g.IrregBlockDouble()
            blk.t = np.hstack([b.t for b in blocks_in])
            for f in fields:
                blk.data[f] = np.hstack([b.data[f.split('.')[-1]] for b in blocks_in])
            blocks_out.append((group_name, blk))
        if raw:
            return blocks_out
        # Reformat for sisock.
        data = {}
        timelines = {}
        for group_name, block in blocks_out:
            timelines[group_name] = {
                't': block.t,
                'finalized_until': block.t[-1],
                'fields': list(block.data.keys()),
            }
            for k,v in block.data.items():
                data[k] = v
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
        """
        unpack = isinstance(fields, str)
        if unpack:
            fields = [fields]
        data, timelines = self.get_data(fields, start, end, min_stride, raw, short_match)
        output = {}
        for timeline in timelines.values():
            # Make the array here, so that the same object is returned
            # for all fields in this group.
            _t = np.array(timeline['t'])
            for f in timeline['fields']:
                output[f] = (_t, np.array(data[f]))
        output = [output[f] for f in fields]
        if unpack:
            output = output[0]
        return output


class _HKProvider:
    def __init__(self, prov_id, prefix):
        self.prov_id = prov_id
        self.prefix = prefix
        self.blocks = []

    @classmethod
    def from_g3(cls, element):
        prov_id = element['prov_id'].value
        prefix = element['description'].value
        return cls(prov_id, prefix)


class HKArchiveScanner:
    """Consume SO Housekeeping streams and create a record of what fields
    cover what time ranges.  This can run as a G3Pipeline module, but
    will not be able to record stream indexing information in that
    case.  If it's populated through the process_file method, then
    index information (in the sense of filenames and byte offsets)
    will be stored.

    """
    def __init__(self):
        self.session_id = None
        self.providers = {}
        self.catalog = HKArchive()
        self.field_groups = []
        self.frame_info = []
        self.counter = -1

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
            
        if f.type == core.G3FrameType.EndProcessing:
            return [f]

        if f.type != core.G3FrameType.Housekeeping:
            return f

        if f['hkagg_type'] == so3g.HKFrameType.session:
            session_id = f['session_id']
            if self.session_id is not None:
                if self.session_id != session_id:
                    self.flush()
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
            representatives = [block['fields'][0] for block in prov.blocks]
            for b in f['blocks']:
                fields = b.data.keys()
                if len(b.t) == 0 or len(fields) == 0:
                    continue
                for block_index,rep in enumerate(representatives):
                    if rep in fields:
                        break
                else:
                    block_index = len(prov.blocks)
                    prov.blocks.append({'fields': fields,
                                        'start': b.t[0],
                                        'index_info': []})
                # To ensure that the last sample is actually included
                # in the semi-open intervals we use to track frames,
                # the "end" time has to be after the final sample.
                prov.blocks[block_index]['end'] = b.t[-1] + SPAN_BUFFER_SECONDS
                prov.blocks[block_index]['index_info'].append(index_info)
                
        else:
            core.log_warn('Weird hkagg_type: %i' % f['hkagg_type'],
                          unit='HKScanner')
        return [f]

    def flush(self, provs=None):
        """Convert any cached provider information into _FieldGroup
        information.  Delete the provider information.  This will be
        called automatically during frame processing if a provider
        session ends.  Once frame processing has completed, this
        fnuction should be called, with no arguments, to flush any
        remaining provider sessions.

        Args:
            provs (list or None): list of provider identifiers to
              flush.  If None, flush all.

        Returns:
            None
        """
        if provs is None:
            provs = list(self.providers.keys())
        for p in provs:
            prov = self.providers.pop(p)
            blocks = prov.blocks
            for info in blocks:
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

    def process_file(self, filename):
        """Process the file specified by ``filename`` using a G3IndexedReader.
        Each frame from the file is passed to self.Process, with the
        optional index_info argument set to a dictionary containing
        the filename and byte_offset of the frame.

        """
        reader = so3g.G3IndexedReader(filename)
        while True:
            info = {'filename': filename,
                    'byte_offset': reader.Tell()}
            frames = reader.Process(None)
            assert len(frames) <= 1
            if len(frames) == 0:
                break
            self(frames[0], info)


class _FieldGroup:
    """Track a set of readout fields that share a timeline.  Attributes
    are:

    - fields (list of str): the field names.
    - cover (IntervalsDouble): time range over which these fields have data.
    - index_info (dict): description of how to find the actual data
      (perhaps a filename and byte_offset?).

    """
    def __init__(self, prefix, fields, start, end, index_info):
        self.prefix = prefix
        self.fields = list(fields)
        self.cover = so3g.IntervalsDouble().add_interval(start, end)
        self.index_info = index_info

if __name__ == '__main__':
    import sys

    hkcs = HKArchiveScanner()
    for filename in sys.argv[1:]:
        hkcs.process_file(filename)
    cat = hkcs.finalize()
    # Get list of fields, timelines, spanning all times:
    fields, timelines = cat.get_fields()
    # Show
    for field_name in sorted(fields):
        group_name = fields[field_name]['timeline']
        print(field_name, timelines[group_name]['interval'])
    full_name = list(fields.keys())[0]
    print('Pretending interest in:', full_name)

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

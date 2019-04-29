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

class HKArchive:
    """Contains information necessary to determine what data fields are
    present in a data archive at what times.  It also knows how to
    group fields that share a commong timeline.
    """
    channel_groups = None
    def __init__(self, channel_groups=None):
        if channel_groups is not None:
            self.channel_groups = list(channel_groups)

    def _get_groups(self, fields=None, start=None, end=None):
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

        Returns:
            List of (group_name, group_fields, cgs).  The
            ``group_name`` is internally generated... probably
            something like "group0".  ``group_fields`` is a list of
            fields belonging to this group.  ``cgs`` is a list of
            _ChannelGroup objects relevant to the fields in this
            group.
        """
        span = so3g.IntervalsDouble()
        if start is None:
            start = span.domain[0]
        if end is None:
            end = span.domain[1]
        span.add_interval(start, end)
        field_map = {}
        for cg in self.channel_groups:
            both = span.intersect(cg.cover)
            if len(both.array()) > 0:
                for f in cg.channels:
                    if fields is not None and f not in fields:
                        continue
                    if not f in field_map:
                        field_map[f] = [cg]
                    else:
                        field_map[f].append(cg)

        # Sort each list of channel_groups, for subsequent comparison.
        [f.sort() for f in field_map.values()]

        # Now group together fields if they have identical
        # channel_group lists (because when they do, they can share a
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
        """Get list of fields that might have a sample in the time interval [start,end).

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
                 raw=False):
        """Load data from specified field(s) between specified times.

        Arguments ``field``, ``start``, ``end`` are as described in
        _get_groups.

        Returns:
            Pair of dictionaries, (data, timelines).  The ``data``
            dictionary is a simple map from channel name to a numpy
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
        grouped = self._get_groups(field, start, end)
        handles = {}  # filename -> G3IndexedReader map.
        blocks_out = []
        for group_name, channels, chgrps in grouped:
            blocks_in = []
            for cg in chgrps:
                for r in cg.index_info:
                    f,o = r['filename'], r['byte_offset']
                    if not f in handles:
                        handles[f] = so3g.G3IndexedReader(f)
                    handles[f].Seek(o)
                    f = handles[f].Process(None)
                    assert(len(f) == 1)
                    # Find the right block.
                    for blk in f[0]['blocks']:
                        if channels[0] in blk.data.keys():
                            blocks_in.append(blk)
                            break
            # Create a new Block for this group.
            blk = so3g.IrregBlockDouble()
            blk.t = np.hstack([b.t for b in blocks_in])
            for c in channels:
                blk.data[c] = np.hstack([b.data[c] for b in blocks_in])
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
        self.channel_groups = []
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
            # If a provider has disappeared, flush its information into a ChannelGroup.
            now_prov_id = [p['prov_id'].value for p in f['providers']]
            to_flush = [p for p in self.providers
                        if not p in now_prov_id]
            for p in to_flush:
                self.flush([p])
                    
            # New providers?
            for p in now_prov_id:
                blocks = self.providers.get(p)
                if blocks is None:
                    self.providers[p] = []

        elif f['hkagg_type'] == so3g.HKFrameType.data:
            # Data frame -- merge info for this provider.
            blocks = self.providers[f['prov_id']]
            representatives = [block['channels'][0] for block in blocks]
            for b in f['blocks']:
                channels = b.data.keys()
                if len(b.t) == 0 or len(channels) == 0:
                    continue
                for block_index,rep in enumerate(representatives):
                    if rep in channels:
                        break
                else:
                    block_index = len(blocks)
                    blocks.append({'channels': channels,
                                   'start': b.t[0],
                                   'index_info': []})
                # To ensure that the last sample is actually included
                # in the semi-open intervals we use to track frames,
                # the "end" time has to be after the final sample.
                blocks[block_index]['end'] = b.t[-1] + SPAN_BUFFER_SECONDS
                blocks[block_index]['index_info'].append(index_info)
                
        else:
            core.log_warn('Weird hkagg_type: %i' % f['hkagg_type'],
                          unit='HKScanner')
        return [f]

    def flush(self, provs=None):
        """Convert any cached provider information into _ChannelGroup
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
            blocks = self.providers.pop(p)
            for info in blocks:
                cg = _ChannelGroup(info['channels'], info['start'],
                                   info['end'], info['index_info'])
                self.channel_groups.append(cg)

    def finalize(self):
        """Finalize the processing by calling flush(), then return an
        HKArchive with all the _ChannelGroup information from this
        scanning session.
        """
        self.flush()
        return HKArchive(self.channel_groups)

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


class _ChannelGroup:
    """Track a set of readout channels that share a timeline.  Attributes
    are:

    - channels (list of str): the field names.
    - cover (IntervalsDouble): time range over which these channels have data.
    - index_info (dict): description of how to find the actual data
      (perhaps a filename and byte_offset?).

    """
    def __init__(self, channels, start, end, index_info):
        self.channels = list(channels)
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
    field_name = list(fields.keys())[0]
    fields, timelines = cat.get_data([field_name])
    x = list(timelines.values())[0]['t']
    y = fields[field_name]
    
    import pylab as pl
    pl.plot(x, y)
    pl.title(field_name)
    pl.show()

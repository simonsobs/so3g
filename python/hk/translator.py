"""Backwards compatibility for older SO HK schemas.

"""

import so3g
import so3g.hk
from spt3g import core


class HKTranslator:
    """Translates SO Housekeeping frames from schema version 0 to version
    1; passes v1 frames through unaffected.  This can be used in a
    G3Pipeline to condition HK streams for processing by v1-compatible
    code.

    Version 1 is not a strict superset of version 0, but the two main
    structural features eliminated in v1 (block names and field
    prefixes) were not really used.

    """
    TARGET_VERSION = 1

    def __init__(self):
        self.stats = {'n_hk': 0,
                      'n_other': 0,
                      'versions': {}}

    def Process(self, f):
        """Translates one frame to the target schema.  Irrelevant frames are
        passed through unmodified.

        Args:
          f: a G3Frame

        Returns:
          A list containing only the translated frame.  G3Pipeline
          compatibility would permit us to return a single frame here,
          instead of a length-1 list.  But we also sometimes call
          Process outside of a G3Pipeline, where a consistent output
          type is desirable.  Returning lists is most
          future-compatible; consumers that want to assume length-1
          should assert it to be true.

        """
        if f.type == core.G3FrameType.EndProcessing:
            core.log_info(str(self.stats))
            return [f]

        if f.type != core.G3FrameType.Housekeeping:
            self.stats['n_other'] += 1
            return [f]

        # It is an HK frame.
        hkagg_version = f.get('hkagg_version', 0)

        self.stats['n_hk'] += 1
        self.stats['versions'][hkagg_version] = self.stats['versions'].get(hkagg_version, 0) + 1

        if hkagg_version == 1:
            # Good enough.
            return [f]

        # Always update the version, even if that's our only change...
        if 'hkagg_version' in f:
            f['hkagg_version_orig'] = hkagg_version
            del f['hkagg_version']
        f['hkagg_version'] = 1

        # No difference in Session/Status for v0 -> v1.
        if f.get('hkagg_type') != so3g.HKFrameType.data:
            return [f]

        # Pop the data blocks out of the frame.
        orig_blocks = f.pop('blocks')
        f['blocks'] = core.G3VectorFrameObject()

        # Now process the data blocks.
        for block in orig_blocks:
            new_block = core.G3TimesampleMap()
            new_block.times = so3g.hk.util.get_g3_time(block.t)
            for k in block.data.keys():
                v = block.data[k]
                new_block[k] = core.G3VectorDouble(v)
            f['blocks'].append(new_block)
        return [f]

    def __call__(self, *args, **kwargs):
        return self.Process(*args, **kwargs)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        usage='This program can be used to convert SO HK Frames to the '
        'latest schema version.')
    parser.add_argument('--output-file', '-o', default='out.g3')
    parser.add_argument('files', nargs='+', help=
                        "SO Housekeeping files to convert.")
    args = parser.parse_args()

    # Run me on a G3File containing a Housekeeping stream.
    core.set_log_level(core.G3LogLevel.LOG_INFO)

    print(f'Streaming to {args.output_file}')
    p = core.G3Pipeline()
    p.Add(core.G3Reader(args.files))
    p.Add(HKTranslator())
    p.Add(core.G3Writer(args.output_file))
    p.Run()

"""Backwards compatibility for older SO HK schemas.

"""

import so3g
import so3g.hk
from spt3g import core


class HKTranslator:
    """Translates SO Housekeeping frames from schema versions {v0, v1} to
    schema version 2.  Passes v2 (or newer) frames through,
    unmodified.  This can be used in a G3Pipeline to condition
    archival HK streams for processing by v2-compatible code.  (Note
    that code that works with the short-lived v1 schema should also
    work on a v2 stream, unless it explicitly rejects based on
    hkagg_version.)

    Version 1/2 are not a strict superset of version 0, but the main
    structural features eliminated in v1 (field prefixes) was not
    really used.

    """

    def __init__(self, target_version=2, future_tolerant=True):
        """Arguments:

          target_version (int): 0, 1, or 2.  Version to which to translate
            the stream.  The code is not able to downgrade a stream.  See
            future_tolerant parameter.

          future_tolerant (bool): Determines the behavior of the
            translator should it encounter a frame with hkagg_version
            higher than target_version.  If future_tolerant is True,
            the frame will be passed through unmodified.  Otherwise, a
            ValueError will be raised.

        """
        self.stats = {'n_hk': 0,
                      'n_other': 0,
                      'versions': {}}
        self.target_version = target_version
        self.future_tolerant = future_tolerant

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
        orig_version = f.get('hkagg_version', 0)

        self.stats['n_hk'] += 1
        self.stats['versions'][orig_version] = self.stats['versions'].get(orig_version, 0) + 1

        if orig_version > self.target_version and not self.future_tolerant:
            raise ValueError(
                ('Translator to v%i encountered v%i, but future_tolerant=False.')
                % (self.TARGET_VERSION, orig_version))

        if orig_version >= self.target_version:
            return [f]

        # Always update the version, even if that's our only change...
        if 'hkagg_version' in f:
            if 'hkagg_version_orig' not in f:
                f['hkagg_version_orig'] = orig_version
            del f['hkagg_version']
        f['hkagg_version'] = self.target_version

        # No difference in Session/Status for v0, v1, v2.
        if f.get('hkagg_type') != so3g.HKFrameType.data:
            return [f]

        if self.target_version == 0:
            return [f]

        if orig_version == 0:
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

        if self.target_version == 1:
            return [f]

        if orig_version <= 1:
            # Add 'block_names'.  Since we don't want to start
            # caching Block Stream information, just compute a good
            # block name based on the alphabetically first field in
            # the block.
            block_names = []
            for block in f['blocks']:
                field_names = list(sorted(block.keys()))
                block_names.append('block_for_%s' % field_names[0])
                assert(len(block_names[-1]) < 256)  # What have you done.
            orig_block_names = []
            f['block_names'] = core.G3VectorString(block_names)

        return [f]

    def __call__(self, *args, **kwargs):
        return self.Process(*args, **kwargs)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        usage='This program can be used to convert SO HK Frames to the '
        'latest schema version.')
    parser.add_argument('--output-file', '-o', default='out.g3')
    parser.add_argument('--target-version', type=int)
    parser.add_argument('files', nargs='+', help=
                        "SO Housekeeping files to convert.")
    args = parser.parse_args()

    # Run me on a G3File containing a Housekeeping stream.
    core.set_log_level(core.G3LogLevel.LOG_INFO)

    translator_args = {}
    if args.target_version is not None:
        translator_args['target_version'] = args.target_version

    print(f'Streaming to {args.output_file}')
    p = core.G3Pipeline()
    p.Add(core.G3Reader(args.files))
    p.Add(HKTranslator(**translator_args))
    p.Add(core.G3Writer(args.output_file))
    p.Run()

import unittest

import so3g
from spt3g import core

import time
import numpy as np

test_file = 'hk_out.g3'

class TestHKSessionHelper(unittest.TestCase):
    """TestCase for so3g.hk, session.py stuff.  The HKSession and
    HKScanner units tests could both live in here, as each is best way
    to run tests on the other.

    """

    def test_00_basic(self):
        """Write a stream of HK frames and scan it for errors."""

        # Write a stream of HK frames.
        # (Inspect the output with 'spt3g-dump hk_out.g3 so3g'.)
        print('Streaming to %s' % test_file)
        w = core.G3Writer(test_file)

        # Create something to help us track the aggregator session.
        hksess = so3g.hk.HKSessionHelper(session_id=None,
                                         hkagg_version=2,
                                         description="Test HK data.")

        # Register a data provider.
        prov_id = hksess.add_provider(
            description='Fake data for the real world.')

        # Start the stream -- write the initial session and status frames.
        w.Process(hksess.session_frame())
        w.Process(hksess.status_frame())

        # Add a bunch of data frames
        t_next = time.time()
        for i in range(10):
            f = hksess.data_frame(prov_id=prov_id, timestamp=t_next)
            hk = core.G3TimesampleMap()
            speed = [1.2, 1.2, 1.3, 1.2, 1.3]
            hk.times = [core.G3Time(_t * core.G3Units.second)
                        for _t in t_next + np.arange(len(speed))]
            hk['position'] = core.G3VectorDouble(np.arange(len(speed)))
            hk['speed'] = core.G3VectorDouble(speed)
            hk['error_bits'] = core.G3VectorInt([10] * len(speed))
            hk['mode_str'] = core.G3VectorString(['ok'] * len(speed))
            t_next += len(hk)
            f['blocks'].append(hk)
            f['block_names'].append('main_block')
            w.Process(f)

        w.Flush()
        del w

        print('Stream closed.\n\n')

        # Now play them back...
        print('Reading back:')
        for f in core.G3File(test_file):
            ht = f.get('hkagg_type')
            if ht == so3g.HKFrameType.session:
                print('Session: %i' % f['session_id'])
            elif ht == so3g.HKFrameType.status:
                print('  Status update: %i providers' % (len(f['providers'])))
            elif ht == so3g.HKFrameType.data:
                print('  Data: %i blocks' % len(f['blocks']))
                for i, block in enumerate(f['blocks']):
                    print('    Block %i' % i)
                    for k, v in block.items():
                        print('    %s' % k, v)

        # Scan and validate.
        print()
        print('Running HKScanner on the test data...')
        scanner = so3g.hk.HKScanner()
        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader(test_file))
        pipe.Add(scanner)
        pipe.Run()

        print('Stats: ', scanner.stats)
        print('Providers: ', scanner.providers)

        self.assertEqual(scanner.stats['concerns']['n_error'], 0)
        self.assertEqual(scanner.stats['concerns']['n_warning'], 0)


if __name__ == '__main__':
    unittest.main()

import unittest
import os

import so3g
from spt3g import core


def write_example_file(filename='hk_out.g3'):
    """Generate an example g3 file with Wiring Frame.

    This is based on other test file writing functions, but is unique in that
    it contains a wiring frame for use in testing Seek/Tell.

    Structure of frames in this file should be:
        - Housekeeping
        - Housekeeping
        - Wiring
        - Housekeeping
        - Housekeeping

    Args:
        filename (str): filename to write data to

    """
    test_file = filename

    # Write a stream of HK frames.
    w = core.G3Writer(test_file)

    # Create something to help us track the aggregator session.
    hksess = so3g.hk.HKSessionHelper(session_id=1234,
                                     description="Test HK data.")

    # Register a data provider.
    prov_id = hksess.add_provider(
        description='Fake data for the real world.')

    # Start the stream -- write the initial session and status frames.
    f = hksess.session_frame()
    w.Process(f)
    f = hksess.status_frame()
    w.Process(f)

    # Write dummy wiring frame
    f = core.G3Frame()
    f.type = core.G3FrameType.Wiring
    w.Process(f)

    # Now make a data frame.
    f = hksess.data_frame(prov_id=prov_id)

    # Add a data block.
    hk = so3g.IrregBlockDouble()
    hk.prefix = 'hwp_'
    hk.data['position'] = [1, 2, 3, 4, 5]
    hk.data['speed'] = [1.2, 1.2, 1.2, 1.2, 1.2]
    hk.t = [0, 1, 2, 3, 4]
    f['blocks'].append(hk)

    # Write two more housekeeping frames.
    w.Process(f)
    w.Process(f)

    del w


class TestG3IndexedReader(unittest.TestCase):
    """TestCase for testing the so3g.G3IndexedReader, which has seek
    capabilities to jump to known frames within a .g3 file.

    """
    def setUp(self):
        self._file = 'test.g3'
        write_example_file(self._file)

    def tearDown(self):
        """Remove the temporary data file we wrote."""
        os.remove(self._file)

    def test_seek(self):
        """Test the Seek/Tell functionality of the G3IndexedReader. We read the
        first four frames, recording the position of the only Wiring frame in the file
        with Tell(). Then we Seek to that location and start reading again, expecting
        the first frame after Seek() to be the wiring frame.

        """
        print("Testing Seek/Tell in G3IndexedReader")
        r = so3g.G3IndexedReader(self._file)
        # Limit the number of Process calls, if we hit the end of the file,
        # then Seek won't work...
        for i in range(4):
            pos = r.Tell()
            f = r.Process(None)[0]
            print("  " + str(f.type))
            if f.type == core.G3FrameType.Wiring:
                w_pos = pos
                print('  Saved wiring frame position: {}'.format(w_pos))

        r.Seek(w_pos)

        # Now that we've seeked, our next frame should be Wiring
        assert r.Process(None)[0].type == core.G3FrameType.Wiring

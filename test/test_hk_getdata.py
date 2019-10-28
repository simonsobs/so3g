import unittest
import os
import numpy as np

import so3g
from so3g.hk import HKArchiveScanner
from so3g.hk import HKTranslator

from spt3g import core

class Seeder(list):
    """Module that can be pre-initialized with frames to feed into a
    Pipeline.  If it's not the first module, it inserts the frames
    prior to passing any others on.

    """
    def Process(self, frames):
        output = [x for x in self]
        self.clear()
        if frames is not None:
            output.extend([f for f in frames if f is not None])
        return output
    def __call__(self, *args, **kw):
        return self.Process(*args, **kw)

def write_example_file(filename='hk_out.g3'):
    """Generate some example HK data and write to file.

    Args:
        filename (str): filename to write data to

    """
    test_file = filename

    # Write a stream of HK frames.
    # (Inspect the output with 'spt3g-dump hk_out.g3 so3g'.)
    seeder = Seeder()
    w = core.G3Pipeline()
    w.Add(seeder)
    w.Add(HKTranslator)
    w.Add(core.G3Writer(test_file))
    

    # Create something to help us track the aggregator session.
    hksess = so3g.hk.HKSessionHelper(session_id=1234,
                                     description="Test HK data.")

    # Register a data provider.
    prov_id = hksess.add_provider(
        description='Fake data for the real world.')

    # Start the stream -- write the initial session and status frames.
    seeder.append(hksess.session_frame())
    seeder.append(hksess.status_frame())

    # Now make a data frame.
    f = hksess.data_frame(prov_id=prov_id)

    # Add a data block.
    hk = so3g.IrregBlockDouble()
    hk.prefix = 'hwp_'
    hk.data['position'] = [1, 2, 3, 4, 5]
    hk.data['speed'] = [1.2, 1.2, 1.2, 1.2, 1.2]
    hk.t = [0, 1, 2, 3, 4]
    f['blocks'].append(hk)

    seeder.append(f)

    w.Run()
    del w


def load_data(filename):
    """Boiled down example of loading the data using an HKArchiveScanner.

    Args:
        filename (str): filename from which to load data

    Returns:
        dict, dict: dictionaries specified in the get_data API

    """
    hkas = HKArchiveScanner()
    hkas.process_file(filename)
    cat = hkas.finalize()
    fields, timelines = cat.get_data(['position'], short_match=True)
    return fields, timelines


class TestGetData(unittest.TestCase):
    """TestCase for testing hk.getdata.py."""
    def setUp(self):
        """Generate some test HK data."""
        self._file = 'test.g3'
        write_example_file(self._file)

    def tearDown(self):
        """Remove the temporary file we made."""
        os.remove(self._file)

    def test_hk_getdata_field_array_type(self):
        """Make sure we return the fields as a numpy array when we get_data."""
        fields, _ = load_data(self._file)
        assert isinstance(fields['position'], np.ndarray)

    def test_hk_getdata_timeline_array_type(self):
        """Make sure we return the timelines as a numpy array when we
        get_data.
        """
        _, timelines = load_data(self._file)
        assert isinstance(timelines['group0']['t'], np.ndarray)


if __name__ == '__main__':
    unittest.main()

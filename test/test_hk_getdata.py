import unittest
import os
import numpy as np

import so3g
from so3g.hk import HKArchiveScanner
from so3g.hk import HKTranslator

from spt3g import core

class Seeder(list):
    """Module that can be pre-initialized with frames to feed into a
    Pipeline.  If it's not the first module, it passes along incoming
    frames and then appends its own.

    """
    def Process(self, frame):
        if frame is not None:
            return [frame]
        output = [x for x in self]
        self.clear()
        return output
    def __call__(self, *args, **kw):
        return self.Process(*args, **kw)

def get_v0_stream():
    """Generate some example HK data, in schema version 0.

    Returns a list of frames constituting a valid version 0 HK stream.

    """
    # Create something to help us track the aggregator session.
    hksess = so3g.hk.HKSessionHelper(session_id=1234,
                                     hkagg_version=0,
                                     description="Test HK data.")

    # Register a data provider.
    prov_id = hksess.add_provider(
        description='Fake data for the real world.')

    # Start the stream -- write the initial session and status frames.
    frames = [
        hksess.session_frame(),
        hksess.status_frame(),
    ]

    # Now make a data frame.
    f = hksess.data_frame(prov_id=prov_id)

    # Add a data block.
    hk = so3g.IrregBlockDouble()
    hk.prefix = 'hwp_'
    hk.data['position'] = [1, 2, 3, 4, 5]
    hk.data['speed'] = [1.2, 1.2, 1.2, 1.2, 1.2]
    hk.t = [0, 1, 2, 3, 4]
    f['blocks'].append(hk)

    frames.append(f)
    return frames

def get_v2_stream():
    """Generate some example HK data, in schema version 2.

    Returns a list of frames constituting a valid version 2 HK stream.

    """
    # Create something to help us track the aggregator session.
    hksess = so3g.hk.HKSessionHelper(session_id=1234,
                                     hkagg_version=2,
                                     description="Test HK data.")

    # Register a data provider.
    prov_id = hksess.add_provider(
        description='Fake data for the real world.')

    # Start the stream -- write the initial session and status frames.
    frames = [
        hksess.session_frame(),
        hksess.status_frame(),
    ]

    # Now make a data frame.
    f = hksess.data_frame(prov_id=prov_id)

    # Add some data blocks.
    hk = core.G3TimesampleMap()
    hk.times = core.G3VectorTime([core.G3Time(i*core.G3Units.seconds) for i in [0, 1, 2, 3, 4]])
    hk['speed'] = core.G3VectorDouble([1.2, 1.2, 1.2, 1.2, 1.2])
    f['blocks'].append(hk)
    f['block_names'].append('group0')

    hk = core.G3TimesampleMap()
    hk.times = core.G3VectorTime([core.G3Time(i*core.G3Units.seconds) for i in [0, 1, 2, 3, 4]])
    hk['position'] = core.G3VectorInt([1, 2, 3, 4, 5])
    hk['mode'] = core.G3VectorString(['going', 'going', 'going', 'going', 'gone/'])
    f['blocks'].append(hk)
    f['block_names'].append('group1')

    frames.append(f)
    return frames


def write_example_file(filename='hk_out.g3', hkagg_version=2):
    """Generate some example HK data and write to file.

    Args:
        filename (str): filename to write data to
        hkagg_version (int): which HK version to write to file

    """
    test_file = filename

    # Write a stream of HK frames.
    # (Inspect the output with 'spt3g-dump hk_out.g3 so3g'.)
    seeder = Seeder()
    w = core.G3Pipeline()
    w.Add(seeder)
    w.Add(HKTranslator(target_version=hkagg_version))
    w.Add(core.G3Writer(test_file))

    if hkagg_version <= 1:
        seeder.extend(get_v0_stream())
    else:
        seeder.extend(get_v2_stream())
    w.Run()
    del w


def load_data(filename, fields=['position', 'speed']):
    """Boiled down example of loading the data using an HKArchiveScanner.

    Args:
        filename (str): filename from which to load data

    Returns:
        dict, dict: dictionaries specified in the get_data API

    """
    hkas = HKArchiveScanner()
    hkas.process_file(filename)
    cat = hkas.finalize()
    fields, timelines = cat.get_data(fields, short_match=True)
    return fields, timelines


class TestGetData(unittest.TestCase):
    """TestCase for testing hk.getdata.py."""
    def setUp(self):
        """Generate some test HK data."""
        self._files = [(0, 'test_0.g3'),
                       (1, 'test_1.g3'),
                       (2, 'test_2.g3')]
        for v, f in self._files:
            write_example_file(f, v)

    def tearDown(self):
        """Remove the temporary file we made."""
        for v, f in self._files:
            os.remove(f)

    def test_hk_getdata_field_array_type(self):
        """Make sure we return the fields as a numpy array when we get_data."""
        for v, f in self._files:
            fields, _ = load_data(f)
            self.assertIsInstance(fields['position'], np.ndarray)

    def test_hk_getdata_timeline_array_type(self):
        """Make sure we return the timelines as a numpy array when we
        get_data.
        """
        for v, f in self._files:
            _, timelines = load_data(f)
            self.assertIsInstance(timelines['group0']['t'], np.ndarray)

    def test_hk_getdata_strings(self):
        """Make sure strings are returned as arrays of unicode."""
        for v, f in self._files:
            if v < 2:
                continue
            fields, _ = load_data(f, ['mode'])
            self.assertIsInstance(fields['mode'], np.ndarray)
            self.assertEqual(fields['mode'].dtype.char, 'U',
                             "String data should unpack to Unicode array.")


if __name__ == '__main__':
    unittest.main()

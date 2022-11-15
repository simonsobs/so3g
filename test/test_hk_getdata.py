import unittest
import os
import numpy as np
import tempfile

import so3g
from so3g.hk import HKArchiveScanner, HKTranslator, HKTree
from so3g.hk.tree import HKRef, HKDeadend

from spt3g import core

# Timestamp near which to place test data.
_BASE_TIME = 1581638400.


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

def get_g3vectortime(dt):
    return core.G3VectorTime((np.asarray(dt) + _BASE_TIME)*core.G3Units.seconds)

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
        description='observatory.test_agent.feeds.test_feed')

    # Start the stream -- write the initial session and status frames.
    frames = [
        hksess.session_frame(),
        hksess.status_frame(),
    ]

    # First data frame ...
    f = hksess.data_frame(prov_id=prov_id)

    # Add some data blocks.
    hk = core.G3TimesampleMap()
    hk.times = get_g3vectortime(range(4))
    hk['speed'] = core.G3VectorDouble([1.2] * 4)
    f['blocks'].append(hk)
    f['block_names'].append('group0')

    hk = core.G3TimesampleMap()
    hk.times = get_g3vectortime(range(5))
    hk['position'] = core.G3VectorInt([1, 2, 3, 4, 5])
    hk['mode'] = core.G3VectorString(['going', 'going', 'going', 'going', 'gone'])
    hk['flag'] = core.G3VectorBool([False, False, True, False, False])
    f['blocks'].append(hk)
    f['block_names'].append('group1')

    # ... save frame.
    frames.append(f)

    # Second data frame ...
    f = hksess.data_frame(prov_id=prov_id)
    hk = core.G3TimesampleMap()
    hk.times = get_g3vectortime(range(3))
    hk['speed'] = core.G3VectorDouble([1.2, 1.2, 1.2])
    f['blocks'].append(hk)
    f['block_names'].append('group0')

    # ... save frame.
    frames.append(f)

    # Third data frame ...
    f = hksess.data_frame(prov_id=prov_id)
    hk = core.G3TimesampleMap()
    hk.times = get_g3vectortime(range(6, 6+6))
    hk['position'] = core.G3VectorInt([6, 7, 8, 9, 10, 11])
    hk['mode'] = core.G3VectorString(['left', 'right', 'left',
                                      'right', 'left', 'halt'])
    hk['flag'] = core.G3VectorBool([False, False, True, False, False, False])
    f['blocks'].append(hk)
    f['block_names'].append('group1')

    # ... save frame.
    frames.append(f)

    return frames, {
        'speed': (np.floating, 7),
        'position': (np.integer, 11),
        'mode': (np.str_, 11),
        'flag': (np.bool_, 11),
    }

def get_v1_stream():
    frames, data = get_v2_stream()
    for f in frames:
        del f['hkagg_version']
        f['hkagg_version'] = 1
        if 'block_names' in f:
            del f['block_names']
    return frames, data

def write_example_file(filename='hk_out.g3', hkagg_version=2):
    """Generate some example HK data and write to file.

    Args:
        filename (str): filename to write data to
        hkagg_version (int): which HK version to write to file

    """
    test_file = os.path.realpath(filename)
    os.makedirs(os.path.split(test_file)[0], exist_ok=True)

    # Write a stream of HK frames.
    # (Inspect the output with 'spt3g-dump hk_out.g3 so3g'.)
    seeder = Seeder()
    w = core.G3Pipeline()
    w.Add(seeder)
    w.Add(HKTranslator(target_version=hkagg_version))
    w.Add(core.G3Writer(test_file))

    if hkagg_version == 1:
        frames, fields = get_v1_stream()
    elif hkagg_version == 2:
        frames, fields = get_v2_stream()

    seeder.extend(frames)
    w.Run()
    del w

    return fields

def load_data(filename, fields=['position', 'speed', 'mode', 'flag'],
              raw=False):
    """Boiled down example of loading the data using an HKArchiveScanner.

    Args:
        filename (str): filename from which to load data

    Returns:
        dict, dict: dictionaries specified in the get_data API

    """
    hkas = HKArchiveScanner()
    hkas.process_file(filename)
    cat = hkas.finalize()
    data = cat.get_data(fields, short_match=True, raw=raw)
    if not raw:
        return data
    else:
        return dict(data)


class TestGetData(unittest.TestCase):
    """TestCase for testing hk.getdata.py."""
    def setUp(self):
        """Generate some test HK data."""
        self._files = [
            [2, 'test_v2.g3', None],
            [1, 'test_v1.g3', None],
        ]
        self._fields = []
        for row in self._files:
            row[2] = write_example_file(row[1], row[0])

    def tearDown(self):
        """Remove the temporary file we made."""
        for v, f, fields in self._files:
            os.remove(f)

    def test_hk_getdata_field_array_type(self):
        """Make sure we return all expected fields, as numpy arrays, with
        expected types and lengths.
        """
        for v, f, fields in self._files:
            data, _ = load_data(f)
            for k, (dtype, n) in fields.items():
                self.assertIsInstance(data[k][0], dtype)
                self.assertEqual(len(data[k]), n)

    def test_hk_getdata_raw(self):
        """Check raw form of get_data."""
        for v, f, fields in self._files:
            data = load_data(f, fields=list(fields.keys()), raw=True)
            for group_name, g3map in data.items():
                for k in g3map.keys():
                    self.assertIsInstance(g3map[k], core.G3FrameObject)
                    self.assertEqual(len(g3map[k]), fields[k][1])

    def test_hk_getdata_timeline_array_type(self):
        """Make sure we return the timelines as a numpy array when we
        get_data.
        """
        for v, f, fields in self._files:
            _, timelines = load_data(f)
            self.assertIsInstance(timelines['group0']['t'], np.ndarray)

    def test_hk_getdata_strings(self):
        """Make sure strings are returned as arrays of unicode."""
        for v, f, fields in self._files:
            if v < 2:
                continue
            data, _ = load_data(f, ['mode'])
            self.assertIsInstance(data['mode'], np.ndarray)
            self.assertEqual(data['mode'].dtype.char, 'U',
                             "String data should unpack to Unicode array.")


class TestHkTree(unittest.TestCase):
    def setUp(self):
        """Generate some test HK data."""
        self._tempdir = tempfile.TemporaryDirectory()
        self._filename = os.path.join(self._tempdir.name,
                                      f'{int(_BASE_TIME/1e5):05d}',
                                      f'{int(_BASE_TIME):05d}.g3')
        self._fields = write_example_file(self._filename)

    def tearDown(self):
        """Remove all temporary files."""
        self._tempdir.cleanup()

    def test_hk_tree(self):
        tree = HKTree(_BASE_TIME-3600, stop=_BASE_TIME+3600,
                      data_dir=self._tempdir.name)
        # Check that attributes resolve as expected
        for ref in [
                tree.test_agent.test_feed,
                tree.test_agent.test_feed.speed]:
            self.assertIsInstance(ref, HKRef)
        for deadend in [
                tree.test_agent.nowhere,
                tree.test_agent.test_feed.nothing]:
            self.assertIsInstance(deadend, HKDeadend)
        # Load a single field
        ref = tree.test_agent.test_feed.mode
        data = ref._load()
        self.assertIsInstance(data, dict)
        self.assertEqual(len(data), 1)
        # Load a bunch of fields
        ref = tree.test_agent
        data2 = ref._load()
        self.assertGreater(len(data2), 1)


if __name__ == '__main__':
    unittest.main()

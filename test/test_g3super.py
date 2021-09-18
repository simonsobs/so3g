import unittest
import time

import numpy as np
import so3g
from spt3g import core


class TestSuperTimestream(unittest.TestCase):

    def test_00_dtypes(self):
        times = core.G3VectorTime((1680000000 + np.arange(1000)) * core.G3Units.s)
        names = ['a', 'b', 'c', 'd', 'e']

        for dtype in ['int32', 'int64', 'float32', 'float64']:
            ts = self._get_ts(4, 100, sigma=0, dtype=dtype)
            ts.encode(1.)
            ts.decode()
            assert(ts.data.dtype is np.dtype(dtype))
            #assert(ts.dtype is np.dtype(dtype))

    def test_01_consistency(self):
        """Test that concordance of (times, names, data) is enforced."""
        names, times, data = self._get_ts(5, 1000, raw=True)

        ts = so3g.G3SuperTimestream()
        ts.times = times
        ts.names = names
        with self.assertRaises(ValueError):
            ts.times = times[:-1]

        with self.assertRaises(ValueError):
            ts.names = names[:-1]

        ts = so3g.G3SuperTimestream()
        with self.assertRaises(ValueError):
            ts.data = data

        ts = so3g.G3SuperTimestream()
        ts.names = names
        with self.assertRaises(ValueError):
            ts.data = data

        ts = so3g.G3SuperTimestream()
        ts.times = times
        with self.assertRaises(ValueError):
            ts.data = data

        ts = so3g.G3SuperTimestream()
        ts.names = names
        ts.times = times
        with self.assertRaises(ValueError):
            ts.data = data[1:]
        with self.assertRaises(ValueError):
            ts.data = data[:,:-1]

    def test_10_encode_int(self):
        ts = self._get_ts(200, 10000)
        d1 = ts.data.copy()
        ts.encode(1)
        np.testing.assert_array_equal(d1, ts.data)

        # Test with a few offsets...
        for offset in [2**25, 2**26 / 3., -1.78 * 2**27]:
            ts = self._get_ts(200, 10000)
            ts.data += int(offset)
            d1 = ts.data.copy()
            ts.encode(1)
            np.testing.assert_array_equal(d1, ts.data)

    def test_20_encode_float(self):
        for dtype in ['float32', 'float64']:
            ts = self._get_ts(200, 10000, sigma=5., dtype=dtype)
            decimals = 2
            precision = 10**(-decimals)
            ts.data[:] = np.round(ts.data, decimals)
            d1 = ts.data.copy()
            ts.encode(10**(-decimals))
            np.testing.assert_allclose(d1, ts.data, atol=precision*1e-3)

    def test_30_serialize(self):
        test_file = 'test_g3super.g3'
        ts = self._get_ts(200, 10000)

        w = core.G3Writer(test_file)
        f = core.G3Frame()
        f['a'] = ts
        w.Process(f)
        del w
        r = core.G3Reader(test_file)
        ts2 = r.Process(None)[0]['a']
        self._check_equal(ts, ts2)

        offsets = [0, 2**25, 2**26 / 3., -1.78 * 2**27]
        records = []
        w = core.G3Writer(test_file)
        for offset in offsets:
            f = core.G3Frame()
            ts.data += int(offset)
            records.append(ts.data.copy())
            f['a'] = ts
            w.Process(f)
        del w
        r = core.G3Reader(test_file)
        for offset, record in zip(offsets, records):
            ts2 = r.Process(None)[0]['a']
            np.testing.assert_array_equal(record, ts2.data)

    # Support functions

    def _get_ts(self, nchans, ntimes, sigma=256, dtype='int32', raw=False):
        names = ['x%i' % i for i in range(nchans)]
        times = core.G3VectorTime(
            (1680000000 + np.arange(ntimes) * .005) * core.G3Units.s)
        data = (np.random.normal(size=(len(names), len(times))) * sigma).astype(dtype)
        if raw:
            return names, times, data
        ts = so3g.G3SuperTimestream()
        ts.names = names
        ts.times = times
        ts.data = data
        return ts

    def _check_equal(self, ts1, ts2):
        np.testing.assert_array_equal(ts1.data, ts2.data)
        np.testing.assert_array_equal(np.array(ts1.times), np.array(ts2.times))
        np.testing.assert_array_equal(np.array(ts1.names), np.array(ts2.names))


def offline_test_memory_leak(MB_per_second=100, encode=True, decode=True):
    """Memory leak loop ... not intended for automated testing!"""
    helper = TestSuperTimestream()
    tick_time = .5
    next_tick = time.time()
    while True:
        d = time.time() - next_tick
        if d < 0:
            time.sleep(d)
        next_tick += tick_time
        print(' ... tick.')
        ts = helper._get_ts(100, int(10000 * MB_per_second * tick_time / 4))
        if encode:
            ts.encode(1)
            if decode:
                ts.decode()

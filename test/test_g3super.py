import unittest

import numpy as np
import so3g
from spt3g import core


class TestSuperTimestream(unittest.TestCase):

    def _get_ts(self, nchans, ntimes, sigma=256):
        chans = ['x%i' % i for i in range(nchans)]
        times = core.G3VectorTime(
            (1680000000 + np.arange(ntimes) * .005) * core.G3Units.s)
        data = (np.random.normal(size=(len(chans), len(times))) * sigma).astype('int32')
        ts = so3g.G3SuperTimestream()
        ts.times = times
        ts.names = chans
        ts.data = data
        return ts

    def _check_equal(self, ts1, ts2):
        np.testing.assert_array_equal(ts1.data, ts2.data)
        np.testing.assert_array_equal(np.array(ts1.times), np.array(ts2.times))
        np.testing.assert_array_equal(np.array(ts1.names), np.array(ts2.names))

    def test_00_basic(self):
        times = core.G3VectorTime((1680000000 + np.arange(10000)) * core.G3Units.s)
        chans = ['a', 'b', 'c', 'd', 'e']
        data = np.zeros((len(chans), len(times)), dtype='int32')

        ts = so3g.G3SuperTimestream()
        ts.times = times
        ts.names = chans
        ts.data = data

        with self.assertRaises(ValueError):
            ts.times = times[:-1]

        with self.assertRaises(ValueError):
            ts.names = chans[:-1]

        ts = so3g.G3SuperTimestream()
        with self.assertRaises(ValueError):
            ts.data = data

        ts = so3g.G3SuperTimestream()
        ts.names = chans
        with self.assertRaises(ValueError):
            ts.data = data

        ts = so3g.G3SuperTimestream()
        ts.times = times
        with self.assertRaises(ValueError):
            ts.data = data

        ts = so3g.G3SuperTimestream()
        ts.names = chans
        ts.times = times
        with self.assertRaises(ValueError):
            ts.data = data[1:]
        with self.assertRaises(ValueError):
            ts.data = data[:,:-1]

    def test_10_encode(self):
        ts = self._get_ts(200, 10000)
        d1 = ts.data.copy()
        ts.encode(1)
        np.testing.assert_array_equal(d1, ts.data)

    def test_20_serialize(self):
        ts = self._get_ts(200, 10000)
        ts.encode(1)
        test_file = 'test_g3super.g3'
        w = core.G3Writer(test_file)
        f = core.G3Frame()
        f['a'] = ts
        w.Process(f)
        del w
        r = core.G3Reader(test_file)
        ts2 = r.Process(None)[0]['a']
        self._check_equal(ts, ts2)

import unittest
import time

import numpy as np
import so3g
from spt3g import core


INT_DTYPES = ['int32', 'int64']
FLOAT_DTYPES = ['float32', 'float64']
ALL_DTYPES = INT_DTYPES + FLOAT_DTYPES


class TestSuperTimestream(unittest.TestCase):

    def test_00_dtypes(self):
        for dtype in ALL_DTYPES:
            ts = self._get_ts(4, 100, sigma=0, dtype=dtype)
            ts.encode()
            ts.decode()
            assert(ts.data.dtype is np.dtype(dtype))

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
        for dtype in INT_DTYPES:
            err_msg = f'Failure during test of dtype={dtype}'
            ts = self._get_ts(10, 980, dtype=dtype)
            d1 = ts.data.copy()
            ts.encode()
            np.testing.assert_array_equal(d1, ts.data, err_msg=err_msg)

            # Test with a few offsets...
            for offset in [2**25, 2**26 / 3., -1.78 * 2**27]:
                ts = self._get_ts(10, 980)
                ts.data += int(offset)
                d1 = ts.data.copy()
                ts.encode()
                np.testing.assert_array_equal(d1, ts.data, err_msg=err_msg)

    def test_20_encode_float(self):
        for dtype in FLOAT_DTYPES:
            err_msg = f'Failure during test of dtype={dtype}'
            ts = self._get_ts(9, 1290, sigma=5., dtype=dtype)
            decimals = 2
            precision=10**-decimals
            ts.data[:] = np.round(ts.data, decimals)
            d1 = ts.data.copy()
            ts.options(precision=precision)
            ts.encode()
            np.testing.assert_allclose(d1, ts.data, atol=precision*1e-3,
                                       err_msg=err_msg)

    def test_40_encoding_serialized(self):
        test_file = 'test_g3super.g3'
        offsets = {
            'int32': [0, 2**25, 2**26 / 3., -1.78 * 2**27],
            'int64': [0, 2**25, 2**26 / 3., -1.78 * 2**27],
            'float32': [0],
            'float64': [0, 2**25, 2**26 / 3., -1.78 * 2**27, 1.8*2**35],
        }
        decimals = 2
        precision=10**-decimals

        w = core.G3Writer(test_file)
        records = []
        for dtype in ALL_DTYPES:
            for offset in offsets[dtype]:
                f = core.G3Frame()
                ts = self._get_ts(4, 100, sigma=100, dtype=dtype)
                ts.data += int(offset)
                if dtype in FLOAT_DTYPES:
                    ts.options(precision=precision)
                    ts.data[:] = ts.data.round(decimals)
                records.append(ts.data.copy())
                f['a'] = ts
                w.Process(f)
        del w
        # readback
        r = core.G3Reader(test_file)
        for dtype in ALL_DTYPES:
            for offset in offsets[dtype]:
                err_msg = f'Failed for dtype={dtype}, offset={offset}'
                ts2 = r.Process(None)[0]['a']
                record = records.pop(0)
                np.testing.assert_allclose(record, ts2.data,
                                           atol=precision*1e-3, err_msg=err_msg)

    def test_50_compression(self):
        test_file = 'test_g3super.g3'

        # Entropy?
        sigma_bits = 8
        sigma = 2**sigma_bits
        _get_ts = lambda dtype: self._get_ts(100, 10000, sigma=sigma, dtype=dtype, seed=12345)

        w = core.G3Writer(test_file)

        sizes = {d: [] for d in ALL_DTYPES}
        for dtype in ALL_DTYPES:
            # No compression
            f = core.G3Frame()
            ts = _get_ts(dtype)
            sizes[dtype].append(ts.data.nbytes)
            if dtype in FLOAT_DTYPES:
                ts.options(precision=1.0)
            ts.options(data_algo=0, times_algo=0)
            f['ts_%s' % dtype] = ts
            w.Process(f)

            # Yes compression
            f = core.G3Frame()
            ts = _get_ts(dtype)
            if dtype in FLOAT_DTYPES:
                ts.options(precision=1.0)
            ts.options(times_algo=0)
            f['ts_%s' % dtype] = ts
            w.Process(f)
        del w

        # Readback
        r = so3g.G3IndexedReader(test_file)
        last = 0
        for dtype in ALL_DTYPES:
            for i in range(2):
                r.Process(None)[0]
                here = r.Tell()
                sizes[dtype].append(here - last)
                last = here

        # Process the results...
        for dtype in ALL_DTYPES:
            err_msg = f'Failed for dtype={dtype}'
            n, s_uncomp, s_comp = sizes[dtype]
            comp_ratio = 1. - (s_uncomp - s_comp) / n
            # But really what matters is the bits-per-word, compressed.
            bits_per_word = comp_ratio * 8 * np.dtype(dtype).itemsize
            #print(dtype, bits_per_word / sigma_bits)
            # I think the theoretical limit is 1.3 or so...
            self.assertLess(bits_per_word, sigma_bits * 1.4, err_msg)

    # Support functions

    def _get_ts(self, nchans, ntimes, sigma=256, dtype='int32', raw=False, seed=None):
        if seed is not None:
            np.random.seed(seed)
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

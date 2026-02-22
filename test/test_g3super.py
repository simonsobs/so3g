import unittest
import collections
import os
import time

import numpy as np
import so3g
from spt3g import core


INT_DTYPES = ['int32', 'int64']
FLOAT_DTYPES = ['float32', 'float64']
ALL_DTYPES = INT_DTYPES + FLOAT_DTYPES


class TestSuperTimestream(unittest.TestCase):

    def test_00_dtypes(self):
        """Test that dtypes supported dtypes are managed properly."""
        return
        for dtype in ALL_DTYPES:
            np_dtype = np.dtype(dtype)
            ts = self._get_ts(4, 100, sigma=0, dtype=dtype)
            assert(ts.data.dtype is np_dtype)
            ts.encode()
            assert(ts.dtype is np_dtype)
            ts.decode()
            assert(ts.data.dtype is np_dtype)
            assert(ts.dtype is np_dtype)

        # Reject big-endian
        with self.assertRaises(ValueError):
            ts = self._get_ts(4, 100, sigma=0, dtype='>i4')

    def test_01_consistency(self):
        """Test that consistency of (times, names, data) shapes is enforced.

        """
        return
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

    def test_02_float_mode(self):
        """Test that rules for entering float mode and setting quanta are enforced."""
        return
        names, times, data_int = self._get_ts(5, 1000, raw=True)
        _, _, data_float = self._get_ts(5, 1000, raw=True, dtype='float32')
        cals = np.ones(len(names))

        def get_base():
            ts = so3g.G3SuperTimestream()
            ts.names = names
            ts.times = times
            return ts

        # Allowed to calibrate int data.
        ts = get_base()
        ts.data = data_int
        ts.calibrate(cals)
        self.assertIs(ts.data.dtype, np.dtype('float32'))

        # Check that calibrate updates cals and data
        x = ts.data.copy()
        ts.calibrate(np.ones(len(cals)) * 3)
        self.assertEqual(ts.quanta[0], 3.)

        # Allowed to set float data after quanta.
        ts = get_base()
        ts.quanta = cals
        ts.data = data_float

        # Not allowed to set float data without first setting quanta.
        ts = get_base()
        with self.assertRaises(ValueError):
            ts.data = data_float

        # Not allowed to set int data after setting quanta.
        ts = get_base()
        ts.quanta = cals
        with self.assertRaises(ValueError):
            ts.data = data_int

        # Not allowed to swap in int data once in float mode
        ts = get_base()
        ts.quanta = cals
        ts.data = data_float
        with self.assertRaises(ValueError):
            ts.data = data_int

        # Prevent passing in wrong number of quanta.
        ts = get_base()
        with self.assertRaises(ValueError):
            ts.quanta = cals[:-1]

    def test_03_idempotency(self):
        """Test that re-encode does nothing, re-decode does nothing."""
        return
        ts = self._get_ts(10, 980, dtype='float32')
        a = ts.data
        self.assertIs(ts.data, a)
        ts.encode()
        # Implicit decode.
        b = ts.data
        self.assertIsNot(b, a)
        # Idempotent decode.
        ts.decode()
        self.assertIs(ts.data, b)

    def test_04_constructors(self):
        # Test int arrays ...
        return
        ts1 = self._get_ts(5, 100, seed=100)
        names, times, data = self._get_ts(5, 100, seed=100, raw=True)
        ts2 = so3g.G3SuperTimestream(names, times)
        ts2.data = data
        self._check_equal(ts1, ts2)
        ts3 = so3g.G3SuperTimestream(names, times, data)
        self._check_equal(ts1, ts3)

        ts4 = so3g.G3SuperTimestream(names, times)
        with self.assertRaises(ValueError):
            ts4.data = data[:,:-1]

        with self.assertRaises(ValueError):
            ts4 = so3g.G3SuperTimestream(names, times, data[:,:-1])

        # Test float arrays ...
        ts1 = self._get_ts(5, 100, seed=100, dtype='float32')
        names, times, data = self._get_ts(5, 100, seed=100, dtype='float32', raw=True)
        quanta = np.ones(len(names))
        ts2 = so3g.G3SuperTimestream(names, times)
        ts2.quanta = quanta
        ts2.data = data
        self._check_equal(ts1, ts2)
        ts3 = so3g.G3SuperTimestream(names, times, data, quanta)
        self._check_equal(ts1, ts3)
        with self.assertRaises(ValueError):
            ts4 = so3g.G3SuperTimestream(names, times, data)

    def test_10_encode_int(self):
        """Test encoding and serialization of integer arrays."""
        return
        for dtype in INT_DTYPES:
            err_msg = f'Failure during test of dtype={dtype}'
            ts = self._get_ts(10, 980, dtype=dtype)
            d1 = ts.data.copy()
            ts.encode()
            np.testing.assert_array_equal(d1, ts.data, err_msg=err_msg)

            # Test with a few offsets...
            for offset in [2**25, 2**26 / 3., -1.78 * 2**27]:
                err_msg1 = f'{err_msg} with offset={offset}'
                ts = self._get_ts(10, 980)
                ts.data += int(offset)
                d1 = ts.data.copy()
                ts.encode()
                np.testing.assert_array_equal(d1, ts.data, err_msg=err_msg1)
                self._readback_compare(ts, err_msg=err_msg1)

    def test_11_fallback(self):
        """Test that the code does not fail to serialize short segments or
        highly random data.

        """
        return
        # Short segments
        for nsamp in range(1, 20):
            ts = self._get_ts(1, nsamp, sigma=0, dtype='int32')
            ts.encode()
            self._readback_compare(ts)
        
        # Random time vector.
        n = 200
        ts = self._get_ts(1, n, sigma=0, dtype='int32')
        ts.times = core.G3VectorTime(
            (np.random.uniform(size=n*8) * 256).astype('uint8').view(dtype='int64'))
        self._readback_compare(ts)

        # Random data array.
        n = 200
        ts = self._get_ts(1, n, sigma=0, dtype='int64')
        ts.data = (np.random.uniform(size=n*8) * 256).astype('uint8').view(dtype='int64').reshape(1,-1)
        self._readback_compare(ts)

        # Small n_det (note 1-10 weren't causing a problem by 11+ were...)
        for n_det in range(1, 20):
            ts = self._get_ts(n_det, 1, sigma=0, dtype='int64')
            ts.data = (np.random.uniform(size=n_det*8) * 256).astype('uint8') \
                .view(dtype='int64').reshape(-1, 1)
            ts.encode()
            self._readback_compare(ts)

    def test_20_encode_float(self):
        return
        for dtype in FLOAT_DTYPES:
            err_msg = f'Failure during test of dtype={dtype}'
            ts = self._get_ts(9, 1290, sigma=5., dtype=dtype)
            precision = .01
            ts.data /= precision
            ts.calibrate([precision] * ts.data.shape[0])
            ts.data[:] = np.round(ts.data)
            d1 = ts.data.copy()
            ts.encode()
            np.testing.assert_allclose(d1, ts.data, atol=precision*1e-3,
                                       err_msg=err_msg)

    def test_30_cpp_interface(self):
        return
        # This is a very basic smoke test.
        ts = so3g.test_g3super(2000, 0, 2000)
        self.assertEqual(ts.data.shape, (3, 2000))
        self.assertTrue(np.all(ts.data[0] == 77.))
        self.assertTrue(np.all(ts.data[1:] == 0.))
        self._readback_compare(ts)
        del ts

        ts = so3g.test_g3super(2000, 100, 1800)
        self.assertEqual(ts.data.shape, (3, 1700))
        self.assertTrue(np.all(ts.data[0] == 77.))
        self.assertTrue(np.all(ts.data[1:] == 0.))
        self._readback_compare(ts)
        del ts
        
    def test_40_encoding_serialized(self):
        return
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
                    ts.data[:] = np.round(ts.data)
                    ts.calibrate([.01] * ts.data.shape[0])
                records.append(ts.data.copy())
                f['a'] = ts
                w.Process(f)

        # show we can write a zero-sample object
        f = core.G3Frame()
        ts = self._get_ts(4, 0)
        f['a0'] = ts
        w.Process(f)

        # flush
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
        return
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
            ts.options(enable=0)
            f['ts_%s' % dtype] = ts
            w.Process(f)

            # Yes compression
            f = core.G3Frame()
            ts = _get_ts(dtype)
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

    def test_60_extract(self):
        """Test selective extraction."""
        return
        def _get_ts():
            ts = self._get_ts(5, 100, seed=100, dtype='float32')
            ts.encode()
            return ts

        ts = _get_ts()
        dest = np.zeros((5, 100), dtype='float32')
        ts.extract(dest)
        np.testing.assert_array_equal(dest, ts.data)

        for dest in [
                np.zeros((4, 100), dtype='float32'),
                np.zeros((100), dtype='float32'),
                np.zeros((5, 100), dtype='int32'),
                ['blech'],
                ]:
            ts = _get_ts()
            with self.assertRaises(ValueError):
                ts.extract(dest)

        idx = np.array([2, 1, 3])
        dest = np.zeros((len(idx), 100), dtype='float32')
        ts = _get_ts()
        ts.extract(dest, None, idx)
        np.testing.assert_array_equal(dest, ts.data[idx])
        ts = _get_ts()
        ts.extract(dest, src_indices=idx)
        np.testing.assert_array_equal(dest, ts.data[idx])

        ts = _get_ts()
        dest = np.zeros((1, 2, 100), dtype='float32')
        with self.assertRaises(ValueError):
            ts.extract(dest, None, idx)

        dest = np.zeros((len(idx), 200), dtype='float32')[:,::2]
        with self.assertRaises(ValueError):
            ts.extract(dest, None, idx)

        # What if decoded already?  This should fail (but not
        # segfault) -- subject to change.
        ts = _get_ts()
        ts.decode()
        idx = np.array([2, 1, 3])
        dest = np.zeros((len(idx), 100), dtype='float32')
        with self.assertRaises(ValueError):
            ts.extract(dest, None, idx)

        # Sample subset
        ts = _get_ts()
        dest = np.zeros((5, 91-10), dtype='float32')
        ts.extract(dest, start=10, stop=91)

        #Injection test.
        ts = _get_ts()
        src_i = np.array([2, 1, 3])
        dest_i = np.array([0, 4, 2])
        dest = np.zeros((5, 91-10), dtype='float32')
        ts.extract(dest, dest_i, src_i, 10, 91)
        np.testing.assert_array_equal(dest[dest_i], ts.data[src_i,10:91])
        self.assertEqual(True, np.all(dest[[1, 3]] == 0))

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
        if dtype in FLOAT_DTYPES:
            ts.quanta = np.ones(len(names))
        ts.data = data
        return ts

    def _check_equal(self, ts1, ts2, err_msg=''):
        for key in ['data', 'times', 'names']:
            np.testing.assert_array_equal(
                np.asarray(getattr(ts1, key)), np.asarray(getattr(ts2, key)),
                err_msg='Fields .%s not equal (detail: %s)' % (key, err_msg))

    def _readback_compare(self, ts, filename='readback_test.g3', cleanup=True,
                          err_msg='(no detail)'):

        """Cache the data from ts, write ts to a file, read it back from file,
        compare to cached data.

        """
        # Cache
        fake_ts = (collections.namedtuple('pseudo_ts', ['times', 'names', 'data']))(
            np.array(ts.times), np.array(ts.names), ts.data.copy())
        # Write
        f = core.G3Frame()
        f['item'] = ts
        core.G3Writer(filename).Process(f)
        # Read
        ts1 = core.G3Reader(filename).Process(None)[0]['item']
        self._check_equal(fake_ts, ts1, err_msg=err_msg)
        if cleanup:
            os.remove(filename)


def offline_test_memory_leak(MB_per_second=100, encode=True, decode=True, dtype='int32',
                             reuse=False):
    """Memory leak loop ... not intended for automated testing!  Pass in
    dtype='int32+' or 'int64+' to trigger float_mode promotoion.

    """
    return
    ts = None
    promotion = False
    if dtype[-1] == '+':
        promotion = True
        dtype = dtype[:-1]
    helper = TestSuperTimestream()
    tick_time = .5
    next_tick = time.time()
    while True:
        d = next_tick - time.time()
        if d > 0:
            time.sleep(d)
        next_tick += tick_time
        print(' ... tick.')
        if ts is None or not reuse:
            ts = helper._get_ts(100, int(10000 * MB_per_second * tick_time / 4), dtype=dtype)
        else:
            ts.data = ts.data.copy()
        if promotion:
            ts.calibrate([1.] * len(ts.data))
        if encode:
            ts.encode()
            if decode:
                ts.decode()

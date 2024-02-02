"""
Test the Ranges and RangesMatrix classes
"""

import unittest
import numpy as np

from so3g.proj import Ranges, RangesMatrix

class TestRanges(unittest.TestCase):
    def test_ranges(self):
        mask = np.array([True, True, False, True, False, True])
        r = Ranges.from_mask(mask)
        self.assertEqual(r.ranges().shape, (3,2))
        np.testing.assert_equal(r.mask(), mask)

    def test_matrix(self):
        f = RangesMatrix.zeros((100,100))
        t = RangesMatrix.ones((100,100))
        np.testing.assert_equal(f[0].mask(), ~(t[0].mask()))
        np.testing.assert_equal(f[0].mask(), (~t[0]).mask())
        np.testing.assert_equal(f[0].mask(), (~t)[0].mask())

        for shape in [(10, 100), (0, 200), (10, 0), (0, 0)]:
            r = RangesMatrix.zeros(shape)
            self.assertEqual(r.shape, shape)
            r = r.copy()
            self.assertEqual(r.shape, shape)
            mask = np.zeros(shape=shape, dtype=bool)
            r = RangesMatrix.from_mask(mask)
            self.assertEqual(r.shape, shape)

    def test_indexing(self):
        r = RangesMatrix.zeros((100, 10000))
        for i in range(len(r)):
            r[i].add_interval(i, i+1000)
        i = np.array([1, 3, 10])
        r1 = r[i]
        for _r, _i in zip(r1, i):
            self.assertEqual(_r.ranges()[0][0], _i)

        s = np.zeros(len(r), bool)
        s[i] = True
        r2 = r[s]
        for _r, _i in zip(r2, i):
            self.assertEqual(_r.ranges()[0][0], _i)

    def test_broadcast(self):
        r0 = RangesMatrix.zeros((100, 1000))
        self.assertEqual(r0.shape, (100, 1000))
        self.assertEqual(r0[None,:,:].shape, (1, 100, 1000))
        self.assertEqual(r0[:,None,:].shape, (100, 1, 1000))

        # It should not be possible to pad or index beyond the
        # outermost dimension.  Ranges isn't very smart about this,
        # but RangesMatrix can be.
        with self.assertRaises(IndexError):
            r0[:,:,None]
        with self.assertRaises(IndexError):
            r0[:,:,0]

    def test_concat(self):
        r0 = RangesMatrix.zeros((10, 100))
        r1 = RangesMatrix.ones ((10, 200))
        r2 = RangesMatrix.ones ((20, 100))

        rc = RangesMatrix.concatenate([r0, r1], axis=1)
        self.assertEqual(rc.shape, (10, 300))
        self.assertEqual(rc[0].mask().sum(), r1[0].mask().sum())

        rc = RangesMatrix.concatenate([r0, r2], axis=0)
        self.assertEqual(rc.shape, (30, 100))

        # Zero size is special case
        rx = RangesMatrix.zeros((0, 100))
        rc = RangesMatrix.concatenate([r0, rx], axis=0)
        self.assertEqual(rc.shape, r0.shape)
        rc = RangesMatrix.concatenate([rx, r0], axis=0)
        self.assertEqual(rc.shape, r0.shape)

        rx = RangesMatrix.zeros((10, 0))
        rc = RangesMatrix.concatenate([r0, rx], axis=1)
        self.assertEqual(rc.shape, r0.shape)
        rc = RangesMatrix.concatenate([rx, r0], axis=1)
        self.assertEqual(rc.shape, r0.shape)

        # These should fail due to shape incompat.
        with self.assertRaises(ValueError):
            rc = RangesMatrix.concatenate([r0, r1], axis=0)
        with self.assertRaises(ValueError):
            rc = RangesMatrix.concatenate([r0, r2], axis=1)

    def test_mask(self):
        """Test conversion from/to boolean mask."""
        # Make sure things work for various shapes.
        for shape in [(10, 200),
                      (10, 4, 200),
                      (10, 0, 200),
                      (0, 200),
                      (10, 4, 0),
                      (200)]:
            # Start from a mask.
            m0 = (np.random.uniform(size=shape) > .8)
            rm = RangesMatrix.from_mask(m0)
            m1 = rm.mask()
            self.assertEqual(rm.shape, m0.shape)
            self.assertEqual(np.all(m0 == m1), True)
            print(shape, m0.sum(), rm.shape)

    def test_int_args(self):
        r = Ranges(1000)
        r.add_interval(10, 20)
        r.add_interval(np.int32(30), np.int32(40))
        with self.assertRaises(ValueError):
            r.add_interval(object(), object())
        self.assertEqual(len(r.ranges()), 2)

    def test_close_gaps(self):
        def _get_gapped(size=0):
            r = Ranges(1000)
            r.append_interval_no_check(10, 20)
            r.append_interval_no_check(20+size, 30+size)
            return r
        r = _get_gapped()
        assert(len(r.ranges()) == 2)
        r.close_gaps(0)
        assert(len(r.ranges()) == 1)

        rr = RangesMatrix([_get_gapped(), _get_gapped(10)])
        rr.close_gaps()
        assert(len(rr[0].ranges()) == 1)
        assert(len(rr[1].ranges()) == 2)
        rr.close_gaps(10)
        assert(len(rr[1].ranges()) == 1)

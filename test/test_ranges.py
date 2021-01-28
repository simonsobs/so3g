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
        self.assertCountEqual(r.ranges().shape, (3,2))
        self.assertCountEqual(r.mask(), mask)

    def test_matrix(self):
        f = RangesMatrix.zeros((100,100))
        t = RangesMatrix.ones((100,100))
        self.assertCountEqual(f[0].mask(), ~(t[0].mask()))
        self.assertCountEqual(f[0].mask(), (~t[0]).mask())
        self.assertCountEqual(f[0].mask(), (~t)[0].mask())

    def test_broadcast(self):
        r0 = RangesMatrix.zeros((100, 1000))
        self.assertCountEqual(r0.shape, (100, 1000))
        self.assertCountEqual(r0[None,:,:].shape, (1, 100, 1000))
        self.assertCountEqual(r0[:,None,:].shape, (100, 1, 1000))

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
        self.assertCountEqual(rc.shape, (10, 300))
        self.assertEqual(rc[0].mask().sum(), r1[0].mask().sum())

        rc = RangesMatrix.concatenate([r0, r2], axis=0)
        self.assertCountEqual(rc.shape, (30, 100))

        # These should fail due to shape incompat.
        with self.assertRaises(ValueError):
            rc = RangesMatrix.concatenate([r0, r1], axis=0)
        with self.assertRaises(ValueError):
            rc = RangesMatrix.concatenate([r0, r2], axis=1)

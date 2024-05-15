import unittest

import so3g
import numpy as np


class TestPolyFill(unittest.TestCase):
    """Test the polynomial gap filling."""

    def test_00_basic(self):
        n = 1000
        p0 = [1.0, 0.2, -0.14]
        t = np.arange(1000) * 0.1
        sig0 = np.polyval(p0, t)[None, :]

        ix, nx = 200, 20
        r = so3g.RangesInt32(n)
        r.add_interval(ix, ix + nx)
        s = r.mask()

        BAD_DATA = -100

        for dtype, fillfunc, tolerance in [
            ("float32", so3g.get_gap_fill_poly, 1e-4),
            ("float64", so3g.get_gap_fill_poly64, 1e-10),
        ]:
            sig1 = sig0.astype(dtype)

            # Write the fill into ex array (inplace=False)
            sig = sig1.copy()
            sig[0, s] = BAD_DATA
            ex = np.zeros((s.sum()), dtype=sig.dtype)
            fillfunc([r], sig, 10, 2, False, ex)
            np.testing.assert_allclose(ex, sig1[0, s], rtol=tolerance)
            np.testing.assert_allclose(BAD_DATA, sig[0, s], rtol=tolerance)

            # Write the fill and extract the filled (inplace=True)
            sig = sig1.copy()
            sig[0, s] = BAD_DATA
            ex = np.zeros((s.sum()), dtype=sig.dtype)
            fillfunc([r], sig, 10, 2, True, ex)
            np.testing.assert_allclose(sig1[0, s], sig[0, s], rtol=tolerance)
            np.testing.assert_allclose(BAD_DATA, ex, rtol=tolerance)


class TestBufferWrapper(unittest.TestCase):
    def test_00(self):
        for array_shape, pattern in [
                ((2, 4), (2, 4)),
                ((2, 4), (-1, -1)),
                ((2, 4), (-2,)),
                ((2, 4), (-1, -1, -2,)),
                ((2, 4), (-2, -1, -1)),
                ((2, 4), (-1, -2, -1)),
        ]:
            a = np.zeros(array_shape)
            so3g.test_buffer_wrapper(a, list(pattern))

        for array_shape, pattern in [
                ((2, 4), (-2, -2)),
                ((2, 4), (2,)),
                ((2, 4), (4,)),
                ((2, 4), (2, 3)),
                ((2, 4), (2, 4, -1)),
                ((2, 4), (2, -1, -1)),
                ((2, 4), (-1, -1, -1, -2)),
        ]:
            a = np.zeros(array_shape)
            with self.assertRaises(RuntimeError):
                so3g.test_buffer_wrapper(a, list(pattern))


class TestJumps(unittest.TestCase):
    """
    Test functions used for jump finding.
    """

    def test_00_clean_flag(self):
        a = np.zeros((1, 100), "int32")
        a[:, 10] = 1
        a[:, 47] = 1
        a[:, 0:8] = 1
        a[:, 50:60] = 1
        a[:, 95:] = 1
        b = a.copy(order="C")
        so3g.clean_flag(b, 5)
        flagged = np.where(b)[1]
        # fmt: off
        np.testing.assert_array_equal(
            flagged,
            np.array([0, 1, 2, 3, 4, 5, 6, 7, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 95, 96, 97, 98, 99,]),
        )
        # fmt: on

    def test_01_matched_jump(self):
        a = np.zeros((1, 100), dtype="float32", order="C")
        a[0, 50:] += 1
        b = np.zeros((1, 100), dtype="int32", order="C")
        min_size = np.ones(1, dtype="float32", order="C")*.5
        so3g.matched_jumps(a, b, min_size, 10)
        flagged = np.where(b)[1]
        np.testing.assert_array_equal(flagged, np.arange(46, 53))

    def test_02_scale_jump(self):
        a = np.zeros((1, 100), dtype="float32", order="C")
        a[0, 50:] += 1
        b = np.zeros((1, 100), dtype="float32", order="C")
        so3g.find_quantized_jumps(a, b, np.zeros(1, "float32", "C"), 10, 0.5)
        flagged = np.where(b)[1]
        np.testing.assert_array_equal(flagged, np.arange(50, 60))


if __name__ == "__main__":
    unittest.main()

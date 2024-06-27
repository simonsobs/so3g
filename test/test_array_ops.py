import unittest

import so3g
import numpy as np


class TestPolyFill(unittest.TestCase):
    """Test the polynomial gap filling.

    """
    def test_00_basic(self):
        n = 1000
        p0 = [1., .2, -.14]
        t = np.arange(1000) * .1
        sig0 = np.polyval(p0, t)[None,:]

        ix, nx = 200, 20
        r = so3g.RangesInt32(n)
        r.add_interval(ix, ix+nx)
        s = r.mask()

        BAD_DATA = -100

        for dtype, fillfunc, tolerance in [
                ('float32', so3g.get_gap_fill_poly, 1e-4),
                ('float64', so3g.get_gap_fill_poly64, 1e-10)]:
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

if __name__ == '__main__':
    unittest.main()

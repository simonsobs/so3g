import unittest

import so3g
import numpy as np


class TestFitting(unittest.TestCase):
    """Test fitting noise model."""

    def test_00_fit_noise(self):

        def noise_model(f, p):
            fknee, w, alpha = p[0], p[1], p[2]
            return w * (1 + (fknee / f) ** alpha)

        ndets = 3;
        nsamps = 1024 // 2 + 1 # assume nperseg = 1024 for psd
        dtype = "float32"
        order = "C"

        lowf = 1.
        fwhite = [10., 100.]

        p0 = np.array([10., 2., 0.7]) # fk, w, alpha
        nparams = len(p0)

        tol = 1e-8 # so3g minimization tolerance
        niter = 200*nparams # so3g max iterations
        epsilon = 1e-5 # so3g gradient perturbation epsilon (uncertainty calculation)

        f = np.linspace(0.01, 200., nsamps, dtype=dtype)
        pxx = np.zeros((ndets, nsamps), dtype=dtype, order=order)

        for i in range(ndets):
            pxx[i,:] = noise_model(f, p0)

        so3g_fitout = np.zeros((ndets, nparams),dtype=dtype, order=order)
        so3g_covout = np.zeros((ndets, nparams),dtype=dtype, order=order)

        so3g.fit_noise(f, pxx, so3g_fitout, so3g_covout, lowf, fwhite[0], fwhite[1], tol, niter, epsilon)

        for i in range(ndets):
            residual = np.abs(p0 - so3g_fitout[i]) / so3g_covout[i]
            np.testing.assert_array_less(residual, 1.0)


if __name__ == "__main__":
    unittest.main()
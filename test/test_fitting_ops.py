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
        nsamps = 1024 // 2 + 1 # Assume nperseg = 1024 for psd
        dtype = "float32"
        order = "C"

        lowf = 1.
        fwhite = [10., 100.]

        p0 = np.array([10., 2., 0.7]) # fknee, w, alpha
        nparams = len(p0)

        tol = 1e-8 # so3g minimization tolerance
        niter = 200*nparams # so3g max iterations
        epsilon = 1e-5 # so3g gradient perturbation epsilon

        f = np.linspace(0.01, 200., nsamps, dtype=dtype)
        pxx = np.zeros((ndets, nsamps), dtype=dtype, order=order)

        for i in range(ndets):
            pxx[i,:] = noise_model(f, p0)

        so3g_params = np.zeros((ndets, nparams), dtype=dtype, order=order)
        so3g_sigmas = np.zeros((ndets, nparams), dtype=dtype, order=order)

        so3g.fit_noise(f, pxx, so3g_params, so3g_sigmas, lowf, fwhite[0],
                       fwhite[1], tol, niter, epsilon)

        # Check if fitted parameters deviate from input by more than 1-sigma
        for i in range(ndets):
            residual = np.abs(p0 - so3g_params[i]) / so3g_sigmas[i]
            np.testing.assert_array_less(residual, 1.0)


if __name__ == "__main__":
    unittest.main()
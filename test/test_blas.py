import unittest
import time
import numpy as np

import so3g

class TestBlas(unittest.TestCase):
    def test_smoke(self):
        """Confirm that BLAS is linked in properly working by calling a
        function that uses it.

        """
        ndet, nfreq, nbin = 50, 100000, 2
        nvec = 3
        dtype = 'float32'

        ft = np.ones((ndet, nfreq), dtype)
        bins = np.zeros((nbin, 2), 'int32')

        iD = np.ones((nbin, ndet), dtype) * 2
        iV = np.zeros((nbin, ndet, nvec), dtype)
        s, norm = 1., 1.2

        # Need some bins.
        bins[0] = [0, nfreq//2]
        bins[1] = [nfreq//2, nfreq]

        t0 = time.time()
        so3g.nmat_detvecs_apply(ft, bins, iD, iV, s, norm)
        print('Elapsed: %.6f' % (time.time() - t0))
        self.assertNotEqual(ft[0,0], 1.)
        

if __name__ == '__main__':
    unittest.main()

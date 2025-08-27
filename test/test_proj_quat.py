import unittest

from so3g.proj import quat
from spt3g import core

import time
import numpy as np

DEG = np.pi/180

class TestCoordSys(unittest.TestCase):
    """TestCase for so3g.proj.quat

    """
    convention_pairs = [
        ('iso',    quat.rotation_iso,    quat.decompose_iso),
        ('lonlat', quat.rotation_lonlat, quat.decompose_lonlat),
        ('xieta',  quat.rotation_xieta,  quat.decompose_xieta),
    ]

    def test_00_inversion(self):
        """Test that the basic rotations / decompositions invert properly."""
        test_args = (.1, .2, .3)
        for name, rotation, decompose in self.convention_pairs:
            q = rotation(*test_args)
            check = decompose(q)
            [self.assertAlmostEqual(x, y) for x, y in zip(check, test_args)]

    def test_10_equivalence(self):
        iso_angles    = (.1, .2, .05)
        lonlat_angles = (iso_angles[1], np.pi/2 - iso_angles[0], iso_angles[2])
        q = quat.rotation_lonlat(*lonlat_angles)
        iso_check = quat.decompose_iso(q)
        [self.assertAlmostEqual(x, y) for x, y in zip(iso_angles, iso_check)]

    def test_20_horizon(self):
        """Test that (xi,eta) are parallel to local (az,el)."""
        az0, el0 = 0, 50*DEG
        q_bore = quat.rotation_lonlat(-az0, el0)
        delta = 0.01*DEG
        for xi, eta in [(delta, 0),
                        (0, delta)]:
            q_det = quat.rotation_xieta(xi, eta)
            lon, lat, gamma = quat.decompose_lonlat(q_bore * q_det)
            az, el = -lon, lat
            if xi != 0:
                self.assertGreater(az, az0)
                self.assertAlmostEqual(el, el0)
            if eta != 0:
                self.assertGreater(el, el0)
                self.assertAlmostEqual(az, az0)

    def test_30_xieta(self):
        """Test that xieta decomposition always interprets rotation as
        a gamma angle; this is relevant for "detectors" at xi,eta =
        (0, 0).

        """
        G = 1.
        q = quat.rotation_xieta(0., 0., G)
        _, _, g = quat.decompose_xieta(q)
        assert abs(g - G) < 1e-9


if __name__ == '__main__':
    unittest.main()

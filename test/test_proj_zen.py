import unittest
import itertools
import warnings

from so3g import proj

import numpy as np
import pylab as pl

# Don't require pixell for testing
try:
    from pixell import enmap
    pixell_found = True
except ModuleNotFoundError:
    pixell_found = False

requires_pixell = unittest.skipIf(pixell_found is False, "pixell not found")

DEG = np.pi/180


def _get_polar_orbits(decs=None):
    # Put together a bunch of trajectories that are circles around the
    # north pole, at declinations decs (degrees).
    if decs is None:
        decs = np.arange(-80, 80.1, 10)
    T0 = 1654041600.
    az, el, t = [], [], []
    circle = np.arange(0, 360)
    for _el in decs:
        az.append(circle),
        el.append(circle*0 + _el)
        t.append(circle*0)
    t, az, el = np.hstack(t) + T0, np.hstack(az) * DEG, np.hstack(el) * DEG
    csl = proj.CelestialSightLine.for_lonlat(-az, el)
    fp = proj.FocalPlane.from_xieta([0., .1*DEG], [0, .1*DEG])
    asm = proj.Assembly.attach(csl, fp)
    return len(t), asm


def _get_polar_crossings(angs=None, amp=10., b=0., det_gamma=0.):
    # Put together a bunch of trajectories that cut through the north
    # pole (or approach at impact parameter b in degrees).  The "angs"
    # is a list of trajectory heading angles, in deg (e.g. 0 and 90
    # correspond a line of constant xi and a line of constant eta).
    # The "amp" is the half-length of each trajectory, in DEG.
    # The "det_gamma" is a tweak to the detector angle, in DEG.
    if angs is None:
        angs = np.arange(0, 180., 10.)
    xi0 = np.linspace(-amp, amp, 101) * DEG
    eta0 = 0*xi0 + b * DEG
    xi, eta, gam, t = [], [], [], []
    for ang in angs:
        c, s = np.cos(ang*DEG), np.sin(ang*DEG)
        xi.append(xi0 * c - eta0*s)
        eta.append(eta0 * c + xi0*s)
        gam.append(xi0 * 0)
        t.append(xi0 * 0)
    xi, eta, gam, t = [np.hstack(x) for x in [xi, eta, gam, t]]
    csl = proj.CelestialSightLine.for_lonlat(xi, eta)
    csl.Q = proj.quat.rotation_xieta(xi, eta, gam)
    fp = proj.FocalPlane.from_xieta([0.], [0.], [det_gamma * DEG])  #, .1*DEG], [0, .1*DEG])
    asm = proj.Assembly.attach(csl, fp)
    return len(t), asm


# List zenithal projections along with their sky coverage
# (4pi/full-sky or 2pi/half-sky).

ZEN_PROJS = {
    'arc': '4pi',
    'sin': '2pi',
    'tan': '2pi',
    'zea': '4pi',
}


class TestZenithalProjEng(unittest.TestCase):
    """Test the Projectionist and supporting structures.

    """
    @requires_pixell
    def test_00_zenithal_radius(self):
        """Test radius function."""
        for pix, area in ZEN_PROJS.items():
            res = 1 
            geom = enmap.geometry2(pos=(np.pi/2, 0), res=1*DEG, shape=(361, 361), proj=pix)
            if abs(geom[1].wcs.lonpole - 180) > 1e-3:
                warnings.warn(f'enmap.geometry2 returns lonpole != 180. [{pix}]')
                geom[1].wcs.lonpole = 180

            p = proj.Projectionist.for_geom(geom[0], geom[1])
            test_decs = ([80., 20., -20., -70.])

            for test_dec in test_decs:
                if test_dec <= 0 and area == '2pi':
                    continue
                n, asm = _get_polar_orbits(decs=[test_dec])
                sig = np.ones((2, n), 'float32')
                mm = enmap.enmap(p.to_map(sig, asm, comps='TQU'), wcs=geom[1])
                dec, ra = mm.posmap()
                s = (mm[0] != 0) * np.isfinite(dec)
                map_dec = dec[s].mean() / DEG
                assert(abs(map_dec - test_dec) < 1)

    @requires_pixell
    def test_10_zenithal_near_pole(self):
        """Test polarization response of zenithal projections, near
        pole.  Indirectly this makes sure things are too wonky at r=0.

        """
        for pix, area in ZEN_PROJS.items():
            res = 1/3600 * DEG
            geom = enmap.geometry2(pos=(np.pi/2, 0), res=res, shape=(500, 500), proj=pix)
            if abs(geom[1].wcs.lonpole - 180) > 1e-3:
                warnings.warn(f'enmap.geometry2 returns lonpole != 180. [{pix}]')
                geom[1].wcs.lonpole = 180

            p = proj.Projectionist.for_geom(geom[0], geom[1])
            for det_gamma in [0., 37.]:
                for b_impact in [0., .001, .010]:
                    n, asm = _get_polar_crossings(
                        amp=.05, b=b_impact, det_gamma=det_gamma)
                    sig = np.zeros((1, n), 'float32')

                    # Purity check...
                    mm = enmap.zeros((3,) + geom[0], geom[1])
                    ang = 34 * DEG
                    mm[1] = np.cos(2*ang)
                    mm[2] = np.sin(2*ang)
                    p.from_map(mm, asm, sig, comps='TQU')
                    expected = np.cos(2*(ang - det_gamma * DEG))
                    np.testing.assert_allclose(sig, expected, rtol=0., atol=1e-6)


if __name__ == '__main__':
    unittest.main()

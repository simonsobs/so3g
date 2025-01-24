"""
Test the proj sight-line code against an ephemeris.
"""

import unittest

import so3g
import numpy as np

import ephem
import datetime

# Do not require qpoint for testing... for now.  qpoint is needed for
# any tests that use CelestialSightLine.az_el rather than
# .naive_az_el.
HAS_QPOINT = False
try:
    import qpoint
    HAS_QPOINT = True
except:
    pass

DEG = so3g.proj.DEG
ARCSEC = DEG / 3600.

# A time at which the GMST was approximately 0h.
SIDEREAL_MIDNIGHT_IN_LONDON = 1501299200.0

# Useful sites.
SOBS = so3g.proj.SITES['so']
LONDON = so3g.proj.EarthlySite(0., 51.5, 0.)

class TestProjAstrometry(unittest.TestCase):

    def assertLonLatNear(self, lonlat0, lonlat1, tol=1*ARCSEC):
        """Computes the cartesian points described by lonlat0 and lonlat1, and
        raises an assertion error if they are further apart than tol
        (radians).  The inputs lonlat0 and lonlat1 must have shape
        (..., 2), and broadcasting will be used if the left dimensions
        differ.

        """
        # Get view that is at least 2-d...
        lonlat0 = np.asarray(lonlat0)
        lonlat1 = np.asarray(lonlat1)
        ob = np.broadcast(lonlat0[..., 0], lonlat1[..., 0])
        v = np.zeros(ob.shape + (3,))
        for L, sgn in [(lonlat0, +1), (lonlat1, -1)]:
            v[..., 0] += sgn*(np.cos(L[..., 0]) * np.cos(L[..., 1]))
            v[..., 1] += sgn*(np.cos(L[..., 0]) * np.sin(L[..., 1]))
            v[..., 2] += sgn*np.sin(L[..., 0])
        self.assertLessEqual(np.max((v**2).sum(axis=-1)), tol)

    def assertAlmostEqualMod(self, first, second, modulus, delta=1*ARCSEC):
        diff = (second - first + modulus / 2) % modulus - modulus / 2
        return self.assertLessEqual(abs(diff), delta)

    def test_naive(self):
        """
        Test that the "naive" boresight coordinates aren't too bad.
        """
        site = LONDON
        t0 = SIDEREAL_MIDNIGHT_IN_LONDON
        for year in [0, .25, .5, .75, -30, 30, 60]:
            print('  Time offset %.2f years...' % year)
            t = (t0 + year*365.25*86400)
            az, el = 180*DEG, 39*DEG
            csl = so3g.proj.CelestialSightLine.naive_az_el(
                *[np.array([x]) for x in [t, az, el]],
                site=site)
            ra0, dec0, c, s = csl.coords()[0]
            # Bang that against our quat code...
            rax, decx, gammax = so3g.proj.quat.decompose_lonlat(csl.Q[0])
            self.assertLonLatNear(np.transpose([ra0, dec0]),
                                  np.transpose([rax, decx]))

            ra0 = ra0 % (2*np.pi)
            print('   so3g:    %12.4f %12.4f' %
                  (float(ra0)/DEG, float(dec0)/DEG))

            ra1, dec1 = get_pyephem_radec(az, el, t, site)
            print('   pyephem: %12.4f %12.4f' %
                  (float(ra1)/DEG, float(dec1)/DEG))

            dra = (((ra1 - ra0 + np.pi) % (2*np.pi)) - np.pi) * np.cos(dec0)
            ddec = dec1 - dec0
            dist = (dra**2 + ddec**2)**.5
            print('   delta:    %12.6f deg = %7.3f arcsec' % (
                dist/DEG, dist/ARCSEC))
            self.assertLess(dist, 1*DEG)

    @unittest.skipIf(not HAS_QPOINT, "qpoint not found")
    def test_precise(self):
        """Test that the full precise boresight computation agrees well
        enough with pyephem.

        """
        site = SOBS
        weather = site.typical_weather

        for year in [0, .25, .5, .75, -30, 30, 60]:
            print('  Time offset %.2f years...' % year)
            t = (1577836800 + year*365.25*86400) - 86400 * 600
            az, el = 40*DEG, 50*DEG
            csl = so3g.proj.CelestialSightLine.az_el(
                *[np.array([x]) for x in [t, az, el]],
                site=site, weather=weather)
            ra0, dec0, c, s = csl.coords()[0]
            ra0 = ra0 % (2*np.pi)
            print('   so3g:    %12.4f %12.4f' %
                  (float(ra0)/DEG, float(dec0)/DEG))

            ra1, dec1 = get_pyephem_radec(az, el, t, site, weather)
            print('   pyephem: %12.4f %12.4f' %
                  (float(ra1)/DEG, float(dec1)/DEG))

            dra = (((ra1 - ra0 + np.pi) % (2*np.pi)) - np.pi) * np.cos(dec0)
            ddec = dec1 - dec0
            dist = (dra**2 + ddec**2)**.5
            print('   delta:    %12.6f deg = %7.3f arcsec' % (
                dist/DEG, dist/ARCSEC))
            self.assertLess(dist, 1*ARCSEC)

    @unittest.skipIf(not HAS_QPOINT, "qpoint not found")
    def test_weather(self):
        """Test that different weather conditions produce similar results to
        pyephem.

        """
        site = SOBS
        for weather in [so3g.proj.weather_factory('vacuum'),
                        so3g.proj.weather_factory('toco')]:
            print('  Weather: %s' % weather)
            t = 1577836800.
            az, el = 40*DEG, 50*DEG
            csl = so3g.proj.CelestialSightLine.az_el(
                *[np.array([x]) for x in [t, az, el]],
                site=site, weather=weather)
            ra0, dec0, c, s = csl.coords()[0]
            ra0 = ra0 % (2*np.pi)
            print('   so3g:    %12.4f %12.4f' %
                  (float(ra0)/DEG, float(dec0)/DEG))

            ra1, dec1 = get_pyephem_radec(az, el, t, site, weather)
            print('   pyephem: %12.4f %12.4f' %
                  (float(ra1)/DEG, float(dec1)/DEG))

            dra = (((ra1 - ra0 + np.pi) % (2*np.pi)) - np.pi) * np.cos(dec0)
            ddec = dec1 - dec0
            dist = (dra**2 + ddec**2)**.5
            print('   delta:    %12.6f deg = %7.3f arcsec' % (
                dist/DEG, dist/ARCSEC))
            self.assertLess(dist, 1*ARCSEC)

    def test_trig_interp(self):
        """This test verifies that the trig approximation (a linear
        interpolation table) is working properly.  It just calls a C++
        function that returns the largest deviation in radians.

        """
        deviation = so3g.test_trig(2**14, 0)
        print('Max inverse trig deviation is %.3f arcsec' % (deviation / ARCSEC))
        self.assertLess(deviation, 0.2*ARCSEC)

    @unittest.skipIf(not HAS_QPOINT, "qpoint not found")
    def test_horizon(self):
        """This test is not astrometric so much as a test that coordinate
        systems are accurately implemented; e.g. that the focal plane
        won't be upsidedown after applying the boresight position.

        """
        # Set up a situation where the parallactic rotation is 0.  Any
        # site looking South, but at an elevation above the South
        # celestial pole, will work.  For safety use lat=0.
        site = LONDON
        az, el = 180*DEG, 60*DEG
        t = 1577836800.
        weather = so3g.proj.weather_factory('vacuum')

        for method in [
                so3g.proj.CelestialSightLine.naive_az_el,
                so3g.proj.CelestialSightLine.az_el,
        ]:
            print('Testing parallactic angle for %s' % method)
            csl = method(
                *[np.array([x]) for x in [t, az, el]],
                site=site, weather=weather)
            ra0, dec0, c, s = csl.coords()[0]
            ra0 = ra0 % (2*np.pi)
            # The parallactic angle should be zero.
            gamma = np.arctan2(s, c)
            print('   so3g:    %12.4f %12.4f %12.4f' %
                  (float(ra0)/DEG, float(dec0)/DEG, float(gamma)/DEG))
            self.assertLess(abs(gamma), 1.*DEG)

    def test_focalplane(self):
        # Focal plane
        xi  = np.array([1., 0., 0.]) * DEG
        eta = np.array([0., 0., 1.]) * DEG
        gamma = np.array([45, 0, -45]) * DEG
        fp = so3g.proj.FocalPlane.from_xieta(xi, eta, gamma)

        # Make a CelestialSightLine with some known equatorial pointing.
        r = so3g.proj.quat.rotation_lonlat(35*DEG, 80*DEG, 15*DEG)
        csl = so3g.proj.CelestialSightLine()
        csl.Q = so3g.proj.quat.G3VectorQuat([r])

        # Attach the FocalPlane to the SightLine and retrieve coords
        # of each detector.
        #asmb = so3g.proj.Assembly.attach(csl, fp)
        coords = csl.coords(fp)

        # Compare to a more direct computation, done here.
        for i, (roffset, resp) in enumerate(fp.items()):
            print('Check detector %d:' % i)
            r1 = r * roffset
            print("AAA", type(r1), r1)
            lon0, lat0, gamma0 = so3g.proj.quat.decompose_lonlat(r1)
            print('  - inline computation: ',
                  lon0 / DEG, lat0 / DEG, gamma0 / DEG)
            lon1, lat1, cosg, sing = coords[i][0]
            gamma1 = np.arctan2(sing, cosg)
            print('  - proj computation: ',
                  lon1 / DEG, lat1 / DEG, gamma1 / DEG)
            self.assertLonLatNear([lon0, lat0], [lon1, lat1])
            self.assertAlmostEqualMod(gamma0, gamma1, 360*DEG, 0.01*DEG)


def get_pyephem_radec(az, el, t, site, weather=None):
    """Compute a single ra,dec pointing from the az, el (radians) and t
    (unix timestamp).  site must be an EarthlySite, and weather
    (optional) a Weather object.

    """
    esite = ephem.Observer()
    esite.lat = site.lat * DEG
    esite.long = site.lon * DEG
    if weather is not None:
        weather.apply(esite)
    d = datetime.datetime.fromtimestamp(t, tz=datetime.timezone.utc)
    Xt = d.year, d.month, d.day, d.hour, d.minute, d.second+d.microsecond*1e-6
    esite.date = ephem.date(Xt)
    return esite.radec_of(az, el)


if __name__ == '__main__':
    unittest.main()

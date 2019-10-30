"""
Test the proj sight-line code against an ephemeris.
"""

import unittest

import so3g
import numpy as np

import ephem
import datetime

DEG = so3g.proj.DEG
ARCSEC = DEG / 3600.


class TestProjAstrometry(unittest.TestCase):

    def test_naive(self):
        site = so3g.proj.EarthlySite.get_named('act')

        for year in [0, .25, .5, .75, -30, 30, 60]:
            print('  Time offset %.2f years...' % year)
            t = (1577836800 + year*365.25*86400) - 86400 * 600
            az, el = 40*DEG, 50*DEG
            csl = so3g.proj.CelestialSightLine.naive_az_el(
                *[np.array([x]) for x in [t, az, el]],
                site=site)
            ra0, dec0, c, s = csl.coords()[0]
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
            self.assertTrue(dist < 1*DEG)

    def test_precise(self):
        site = so3g.proj.EarthlySite.get_named('act')
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
            self.assertTrue(dist < 1*ARCSEC)

    def test_weather(self):
        site = so3g.proj.EarthlySite.get_named('act')
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
            self.assertTrue(dist < 1*ARCSEC)


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
    d = datetime.datetime.utcfromtimestamp(t)
    Xt = d.year, d.month, d.day, d.hour, d.minute, d.second+d.microsecond*1e-6
    esite.date = ephem.date(Xt)
    return esite.radec_of(az, el)


if __name__ == '__main__':
    unittest.main()

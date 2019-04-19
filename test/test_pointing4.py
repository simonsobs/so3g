"""
Test the approximate horizon -> celestial pointing code.
"""

import so3g
import numpy as np

import ephem
import time, datetime, calendar

DEG = so3g.proj.DEG

site = so3g.proj.EarthlySite.get_named('act')

for year in [0, .25, .5, .75, -30, 30, 60]:
    t = (1577836800 + year*365.25*86400) - 86400 * 600
    print('Time is %i' % t)
    az, el = 40*DEG, 50*DEG
    csl = so3g.proj.CelestialSightLine.naive_az_el(
        *[np.array([x]) for x in [az,el,t]],
        site=site)

    print(' Getting simple pointing.')
    dets = np.array([(1., 0., 0., 0)])
    X, = csl.coords(dets)
    ra0, dec0 = X[0][:2]
    ra0 = ra0 % (2*np.pi)
    print('   %12.4f %12.4f' % (float(ra0)/DEG,float(dec0)/DEG))

    print(' Check on ephem...')

    esite = ephem.Observer()
    esite.lat = site.lat
    esite.long = site.lon
    d = datetime.datetime.utcfromtimestamp(t)
    Xt = d.year, d.month, d.day, d.hour, d.minute, d.second+d.microsecond*1e-6
    esite.date = ephem.date(Xt)
    ra, dec = esite.radec_of(az, el)
    print('   %12.4f %12.4f' % (float(ra)/DEG,float(dec)/DEG))
    dra = (((ra - ra0 + np.pi) % (2*np.pi)) - np.pi) * np.cos(dec0)
    ddec = dec - dec0
    dist = (dra**2 + ddec**2)**.5
    print('Difference on sky is approximately %.3f degrees.' % (dist/DEG))

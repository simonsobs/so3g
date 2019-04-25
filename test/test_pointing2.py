import so3g
import so3g.proj as SP
from spt3g import core
import numpy as np

from pixell import enmap
from test_utils import Timer
import pylab as pl

from astropy.wcs import WCS
wcs = WCS(naxis=2)
naxis = (128,128)
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs.wcs.crval = np.array([0., 0.])
wcs.wcs.crpix = np.array([naxis[0]//2, naxis[1]//2])[::-1]
wcs.wcs.cdelt = np.array([.1, .1])

# A slightly polarized beam.
#wcs = enmap.wcsutils.cea((0,0), (.1, .1), )
beam = enmap.zeros((3,) + naxis, wcs=wcs)
yx = beam.posmap() * 180/np.pi
beam[:] = np.exp(-(yx[0]**2 + yx[1]**2) / 5**2)[None,...]
phi = np.arctan2(yx[0], yx[1])
beam[1] *= 0.04 * np.cos(2*phi)
beam[2] *= 0.04 * np.sin(2*phi)

def get_pixelizor(emap):
    ny, nx = emap.shape[1:]
    dy, dx = emap.wcs.wcs.cdelt
    y0, x0 = emap.wcs.wcs.crval
    iy0, ix0 = emap.wcs.wcs.crpix
    return so3g.Pixelizor2_Flat(ny, nx, dy, dx, y0, x0, iy0, ix0)

pxz = get_pixelizor(beam)

n_det = 100
n_t = 100000
x = (20 * np.arange(n_t) / n_t) % 1. * 15 - 7.5
y = np.arange(n_t) / n_t * 15 - 7.5

# boresight
ptg = np.zeros((n_t, 4))
ptg[...,0] = x
ptg[...,1] = y
ptg[...,2] = 1.
ptg[...,3] = 0.

# dets... evenly distributed through a ~circle.
r = np.arange(n_det)**.5 / n_det**.5 * .5
ophi = 6.28 * np.arange(n_det)**2 / n_det**2
ofs = np.transpose([r*np.cos(ophi), r*np.sin(ophi), r*0, r*0])
# pol dirs are random.
polphi = 6.28 * np.random.uniform(size=n_det)
ofs[:,2] = np.cos(polphi)
ofs[:,3] = np.sin(polphi)

# Project map into a signal.
pe = so3g.ProjectionEngine2(pxz)

sig0 = pe.from_map(beam, ptg, ofs, None, None)

# Add some noise...
sig1 = sig0 + .0001*np.random.normal(0, 1, size=sig0.shape)

# Then back to map.
map1 = pe.to_map(None, ptg, ofs, sig1, None)

# weights...
wmap1 = pe.to_weight_map(None, ptg, ofs, None, None)
wmap1[1,0] = wmap1[0,1]
wmap1[2,0] = wmap1[0,2]
wmap1[2,1] = wmap1[1,2]
# inverted...
iwmap1 = np.linalg.pinv(wmap1.transpose((2,3,0,1))).transpose((2,3,0,1))

map2 = (iwmap1 * map1[None,...]).sum(axis=1)

pl.rcParams.update({'image.cmap': 'gray'})
for comp in [0,1,2]:
    _, sps = pl.subplots(2, 2)
    sps[0,0].imshow(beam[comp])
    sps[0,1].imshow(map1[comp])
    sps[1,0].imshow(map2[comp])
    sps[1,1].imshow(map2[comp] - beam[comp])
    pl.show()

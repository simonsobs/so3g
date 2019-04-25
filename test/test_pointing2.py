import so3g
import numpy as np
import pylab as pl

from pixell import enmap
from astropy.wcs import WCS

from test_utils import Timer

# Create a world coordinate system, with .1 deg pixels.
wcs = WCS(naxis=2)
naxis = (128,128)
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs.wcs.crval = np.array([0., 0.])
wcs.wcs.crpix = np.array([naxis[0]//2, naxis[1]//2])[::-1]
wcs.wcs.cdelt = np.array([.1, .1])

# A slightly polarized beam.
beam = enmap.zeros((3,) + naxis, wcs=wcs)
yx = beam.posmap() * 180/np.pi
beam[:] = np.exp(-(yx[0]**2 + yx[1]**2) / 5**2)[None,...]
phi = np.arctan2(yx[0], yx[1])
beam[1] *= 0.04 * np.cos(2*phi)
beam[2] *= 0.04 * np.sin(2*phi)

def get_pixelizor(emap):
    # Returns a Pixelizor suitable for use with emap's wcs.
    ny, nx = emap.shape[1:]
    dy, dx = emap.wcs.wcs.cdelt
    y0, x0 = emap.wcs.wcs.crval
    iy0, ix0 = emap.wcs.wcs.crpix
    return so3g.Pixelizor2_Flat(ny, nx, dy, dx, y0, x0, iy0, ix0)

# Use this pixelizor for projections.
pxz = get_pixelizor(beam)

# Our dummy observation.
n_det = 200
n_t = 100000

# Boresight motion.
x = (20 * np.arange(n_t) / n_t) % 1. * 15 - 7.5
y = np.arange(n_t) / n_t * 15 - 7.5

# At each time step, boresight is (x, y, cos(phi), sin(phi))
ptg = np.zeros((n_t, 4))
ptg[...,0] = x
ptg[...,1] = y
ptg[...,2] = 1.
ptg[...,3] = 0.

# Detector offsets... make a diskular bundle of these.
r = np.arange(n_det)**.5 / n_det**.5 * .5
ophi = 6.28 * np.arange(n_det)**2 / n_det**2
ofs = np.transpose([r*np.cos(ophi), r*np.sin(ophi), r*0, r*0])

# Use random polarization directions.
polphi = 6.28 * np.random.uniform(size=n_det)
ofs[:,2] = np.cos(polphi)
ofs[:,3] = np.sin(polphi)

# Actually, make them col-located orthogonal pairs.
ofs[1::2] = ofs[0::2]
ofs[1::2,2] = -ofs[0::2,3]
ofs[1::2,3] =  ofs[0::2,2]

#
# Projection tests.
#

# This is a TQU optimized projector.
pe = so3g.ProjectionEngine2(pxz)

# Project the map into time-domain.
sig0 = pe.from_map(beam, ptg, ofs, None, None)

# Add some noise...
sig1 = sig0 + .0001*np.random.normal(0, 1, size=sig0.shape)

# Then back to map.
map1 = pe.to_map(None, ptg, ofs, sig1, None)

# Get the weight map (matrix).
wmap1 = pe.to_weight_map(None, ptg, ofs, None, None)
wmap1[1,0] = wmap1[0,1]  # fill in unpopulated entries...
wmap1[2,0] = wmap1[0,2]
wmap1[2,1] = wmap1[1,2]

# The inverse weight map.
iwmap1 = np.linalg.pinv(wmap1.transpose((2,3,0,1))).transpose((2,3,0,1))

# The mapped solution.
map2 = (iwmap1 * map1[None,...]).sum(axis=1)

# Show binned map, solved map, residual.
pl.rcParams.update({'image.cmap': 'gray'})
for comp in [0,1,2]:
    _, sps = pl.subplots(2, 2)
    sps[0,0].imshow(beam[comp])
    sps[0,1].imshow(map1[comp])
    sps[1,0].imshow(map2[comp])
    sps[1,1].imshow(map2[comp] - beam[comp])
    pl.show()


#
# Do it again, but with "pair differencing".
#

sigd = (sig1[:,0::2,:] - sig1[:,1::2,:]) / 2
ofsd = (ofs[::2,...])

# Use the QU projector.
pe = so3g.ProjectionEngine1(pxz)

# Bin the map again...
map1d = pe.to_map(None, ptg, ofsd, sigd, None)

# Bin the weights again...
wmap1d = pe.to_weight_map(None, ptg, ofsd, None, None)
wmap1d[1,0] = wmap1d[0,1] # fill in unpopulated entries again...

# Solved...
iwmap1d = np.linalg.pinv(wmap1d.transpose((2,3,0,1))).transpose((2,3,0,1))
map2d = (iwmap1d * map1d[None,...]).sum(axis=1)

# Display.
for comp in [0,1]:
    _, sps = pl.subplots(2, 2)
    sps[0,0].imshow(beam[comp+1])
    sps[0,1].imshow(map1d[comp])
    sps[1,0].imshow(map2d[comp])
    sps[1,1].imshow(map2d[comp] - beam[comp+1])
    pl.show()

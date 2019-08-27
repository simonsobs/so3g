import so3g
import numpy as np
import pylab as pl

pl.rcParams['image.origin'] = 'lower'

from pixell import enmap
from astropy.wcs import WCS

import test_utils
from test_utils import Timer, proj_dict

# Command line option(s)...
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--system', '-s', choices=proj_dict.keys(),
                    default=list(proj_dict.keys())[0])
args = parser.parse_args()
system = args.system
print('Using system: %s' % system)

# Create a world coordinate system, with .1 deg pixels.
wcs = WCS(naxis=2)
naxis = (128,128)
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs.wcs.crval = np.array([0., 0.])
wcs.wcs.crpix = np.array([naxis[0]//2, naxis[1]//2])[::-1]
wcs.wcs.cdelt = np.array([-.1, .1])

# A slightly polarized beam.
beam = enmap.zeros((3,) + naxis, wcs=wcs)
yx = beam.posmap() * 180/np.pi
# Main beam...
beam[:] = np.exp(-(yx[0]**2 + yx[1]**2) / 2**2)[None,...]
# feature, up and to the right.
beam[:] += np.exp(-((yx[0]-3.)**2 + (yx[1]+1.5)**2) / .7**2)[None,...]
phi = np.arctan2(yx[0], yx[1])
beam[1] *= 0.04 * np.cos(2*phi)
beam[2] *= 0.04 * np.sin(2*phi)

def get_pixelizor(emap):
    # Returns a Pixelizor suitable for use with emap's wcs.
    ny, nx = emap.shape[1:]
    dy, dx = emap.wcs.wcs.cdelt * np.pi/180
    iy0, ix0 = emap.wcs.wcs.crpix
    return so3g.Pixelizor2_Flat(ny, nx, dy, dx, iy0, ix0)

# Use this pixelizor for projections.
pxz = get_pixelizor(beam)

# Our dummy observation.
n_det = 200
n_t = 10000#0

# Boresight motion, in degrees, rough projection.
x = (20 * np.arange(n_t) / n_t) % 1. * 15 - 7.5
y = np.arange(n_t) / n_t * 15 - 7.5

# Detector offsets... make a diskular bundle of these.
r = np.arange(n_det)**.5 / n_det**.5 * .5 * (np.pi/180)
ophi = 6.28 * np.arange(n_det)**2 / n_det**2
dx, dy = r*np.cos(ophi), r*np.sin(ophi)
polphi = 6.28 * np.random.uniform(size=n_det)

# ... but not totally random...
polphi[0] = 0.       # The "Q" detector is nice to have for debugging...
polphi[2] = np.pi/4  # ...as is the "U" detector.

# ... and (2n,2n+1) should be co-located orthogonal pairs.
dx[1::2] = dx[0::2]
dy[1::2] = dy[0::2]
polphi[1::2] = polphi[0::2] + np.pi/2

ptg = test_utils.get_boresight_quat(system, x*np.pi/180, y*np.pi/180)
ofs = test_utils.get_offsets_quat(system, dx, dy, polphi)


#
# Projection tests.
#

pe = test_utils.get_proj(system, 'TQU')(pxz)

# Project the map into time-domain.
sig0 = pe.from_map(beam, ptg, ofs, None, None)
sig0 = np.array(sig0)

# Add some noise...
sig1 = sig0 + .0001*np.random.normal(0, 1, size=sig0.shape).astype('float32')

# Show this...
_, sps = pl.subplots(1, 2)
DET = 0
sps[0].plot(sig0[DET] - sig0[DET+1])
sps[1].plot(sig1[DET] - sig1[DET+1])
pl.show()

# Then back to map.
print('Create timestream...', end='\n ... ')
with Timer():
    map1 = pe.to_map(None, ptg, ofs, sig1, None)

# Get the weight map (matrix).
wmap1 = pe.to_weight_map(None, ptg, ofs, None, None)
wmap1[1,0] = wmap1[0,1]  # fill in unpopulated entries...
wmap1[2,0] = wmap1[0,2]
wmap1[2,1] = wmap1[1,2]

# The inverse weight map.
iwmap1 = test_utils.linalg_pinv(wmap1.transpose((2,3,0,1))).transpose((2,3,0,1))

# The mapped solution.
map2 = (iwmap1 * map1[None,...]).sum(axis=1)

# Show binned map, solved map, residual.
pl.rcParams.update({'image.cmap': 'gray'})
mask = wmap1[0,0] > 0
for comp in [0,1,2]:
    _, sps = pl.subplots(2, 2)
    sps[0,0].imshow(beam[comp])
    sps[0,1].imshow(map1[comp])
    sps[1,0].imshow(map2[comp])
    sps[1,1].imshow((map2[comp] - beam[comp])*mask)
    pl.show()


#
# Do it again, but with "pair differencing".
#

sigd = (sig1[0::2,:] - sig1[1::2,:]) / 2
ofsd = (ofs[::2,...])

# Use the QU projector.
pe = test_utils.get_proj(system, 'QU')(pxz)

# Bin the map again...
map1d = pe.to_map(None, ptg, ofsd, sigd, None)

# Bin the weights again...
wmap1d = pe.to_weight_map(None, ptg, ofsd, None, None)
wmap1d[1,0] = wmap1d[0,1] # fill in unpopulated entries again...

# Solved...
iwmap1d = test_utils.linalg_pinv(wmap1d.transpose((2,3,0,1))).transpose((2,3,0,1))
map2d = (iwmap1d * map1d[None,...]).sum(axis=1)

# Display.
mask = wmap1d[0,0] > 0
for comp in [0,1]:
    _, sps = pl.subplots(2, 2)
    sps[0,0].imshow(beam[comp+1])
    sps[0,1].imshow(map1d[comp])
    sps[1,0].imshow(map2d[comp])
    sps[1,1].imshow((map2d[comp] - beam[comp+1])*mask)
    pl.show()

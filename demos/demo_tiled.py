import matplotlib

import so3g
import numpy as np
import pylab as pl

import test_utils
from test_utils import Timer, proj_dict
import os

# Command line option(s)...
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--system', '-s',
                    choices=proj_dict.keys(),
                    default=list(proj_dict.keys())[0])
parser.add_argument('--n-det', '-d',
                    type=int, default=2000)
parser.add_argument('--n-time', '-t',
                    type=int, default=25000)
#parser.add_argument('--tiled', action='store_true')
args = parser.parse_args()
args.tiled = True
system = args.system
print('Using system: %s' % system)

# Map space.
class TiledOccupation:
    def __init__(self, super_shape, tile_shape, maps_present=[]):
        self.super_shape, self.tile_shape = super_shape, tile_shape
        self.maps_present = maps_present
        self.tilage = [int(np.ceil(s/t)) for s, t in zip(self.super_shape, self.tile_shape)]
    def zeros(self, shape):
        if np.ndim(shape) == 0:
            shape = (shape,)
        n = shape[0]
        for s in shape[1:]:
            n *= s
        tiles = pxzt.zeros(n, maps_present)
        for t in tiles:
            if t is not None:
                t.shape = shape + t.shape[1:]
        return tiles
    def new_map(self, shape):
        return [(self.super_shape, self.tile_shape), self.zeros(shape)]
    def __repr__(self):
        return 'TiledOccupation(%s,%s,occupation=%i/%i)' % (
            self.super_shape, self.tile_shape, len(self.maps_present),
            self.tilage[0] * self.tilage[1])
    @staticmethod
    def max_diff(a, b):
        return max([abs(_a-_b).max() for _a, _b in zip(a[1], b[1]) if _a is not None])
        
#super_shape = (598, 502)
#tile_shape = (100, 120)
#tiling = [int(np.ceil(s/t)) for s,t in zip(super_shape, tile_shape)]
tiling = TiledOccupation((598, 502), (100,120))

pxz = so3g.Pixelizor2_Flat(
    tiling.super_shape[0], tiling.super_shape[1],
    0.00005,0.00005,0,0)
pxzt = so3g.Pixelizor2_Flat_Tiled(
    tiling.super_shape[0], tiling.super_shape[1],
    0.00005,0.00005,0,0,
    tiling.tile_shape[0], tiling.tile_shape[1])

# Samples
n_det = args.n_det
n_t = args.n_time

# Boresight
phi = np.arange(n_t) / n_t * 6.28 * 3
r = np.arange(n_t) / n_t * .006
x = r*np.cos(phi) + .004 + .0005 * np.random.uniform(size=len(phi))
y = r*np.sin(phi) + .008 + .0005 * np.random.uniform(size=len(phi))

# offsets
dr = .002
polphi = 6.28 * np.arange(n_det) / n_det
dx, dy = dr * np.cos(polphi), dr * np.sin(polphi)

pe = test_utils.get_proj(system, 'TQU', pxzt, tiled=args.tiled)
ptg = test_utils.get_boresight_quat(system, x, y)
ofs = test_utils.get_offsets_quat(system, dx, dy, polphi)

sig = np.ones((1,n_det,n_t), 'float32') * .5

print('Note the problem size is %i x %i = %.3fM samples.\n' %
      (n_det, n_t, n_det*n_t/1e6))

n_omp = os.getenv('OMP_NUM_THREADS')
if n_omp is None:
    print(' OMP_NUM_THREADS not set -- unknown parallelism.\n')
    n_omp = '?'
else:
    n_omp = int(n_omp)
    print(' OMP_NUM_THREADS=%i.\n' % n_omp)

#
# Projector exercise
#

# 1. The .pixels method should tell us which tiles we will need.

print('Compute coords and pixels and return pixels.', end='\n ... ')
pix = np.empty(sig.shape[1:] + (3,), 'int32')
with Timer() as T:
    spix = pe.pixels(ptg,ofs,pix)

maps_present = []
for s in spix:
    maps_present.extend(list(set(list(s[:,0]))))
    
maps_present = sorted(list(set(maps_present)))
tiling.maps_present = maps_present

# 2. The .zeros method should make tiles of the right sizes.
with Timer() as T:
    tiles1 = pe.to_map(tiling.new_map(3),ptg,ofs,sig[0],None)

print('Compute pixel_ranges (OMP prep)... ', end='\n ... ')
with Timer():
    Ivals = pe.pixel_ranges(ptg, ofs)

if 1:
    coo = np.empty(sig.shape[1:] + (4,), 'double')
    print('Compute and return coordinates only.', end='\n ... ')
    with Timer() as T:
        pe.coords(ptg,ofs[:,:],coo)

    print('And again but with internal creation.', end='\n ... ')
    with Timer() as T:
        coo1 = pe.coords(ptg,ofs,None)
    assert(np.all([(c == _c).all() for c,_c in zip(coo, coo1)]))

    pl.plot(coo[0,:,0],
            coo[0,:,1])
    #pl.show()

    del coo, coo1

    print('Compute coords and pixels and return pixels.', end='\n ... ')
    pix = np.empty(sig.shape[1:] + (3,), 'int32')
    with Timer() as T:
        pe.pixels(ptg,ofs,pix)

    pix[:] = 0
    pix_list = [p for p in pix]  #listify...
    print('And into a list.', end='\n ... ')
    with Timer() as T:
        pe.pixels(ptg,ofs,pix_list)

    print('And with no target array(s).', end='\n ... ')
    with Timer() as T:
        pix3 = pe.pixels(ptg, ofs, None)

    print('Check for failure...', end='\n ... ')
    pix_list[10] = np.empty(sig.shape[-1]*2, 'int32')[::2]
    try:
        pe.pixels(ptg,ofs,pix_list)
        raise Exception('Failed to fail when given invalid data!')
    except RuntimeError as e:
        print('Failure successful.')

    print('Get spin projection factors, too.',
          end='\n ... ')
    with Timer() as T:
        pix2, spin_proj = pe.pointing_matrix(ptg, ofs, None, None)

    del pix, pix_list, pix3, pix2, spin_proj


if 1:
    print('Forward projection (TQU)', end='\n ... ')
    sig_list = [x for x in sig[0]]
    with Timer() as T:
        map1 = pe.to_map(tiling.new_map(3),ptg,ofs,sig_list,None)

if 1:
    print('Forward projection (TQU) with OMP (%s): ' % n_omp, end='\n ... ')
    with Timer() as T:
        map1o = pe.to_map_omp(tiling.new_map(3),ptg,ofs,sig_list,None,Ivals)

if 1:
    print('Reverse projection (TQU)', end='\n ... ')
    sig1 = [x for x in np.zeros((n_det,n_t), 'float32')]
    sig1 = np.random.uniform(size=(n_det,n_t)).astype('float32')
    with Timer() as T:
        #sig1 =
        pe.from_map(map1, ptg, ofs, sig1, None)

if 1:
    print('Forward project weights (TQU)', end='\n ... ')
    map0 = tiling.new_map((3,3))
    with Timer() as T:
        map2 = pe.to_weight_map(map0,ptg,ofs,None,None)

if 1:
    print('Forward project weights (TQU) with OMP (%s): ' % n_omp, end='\n ... ')
    map0 = tiling.new_map((3,3))
    with Timer() as T:
        map2o = pe.to_weight_map_omp(map0,ptg,ofs,None,None,Ivals)

    print('Checking that OMP and non-OMP forward calcs agree: ', end='\n ... ')
    assert tiling.max_diff(map1, map1o) == 0
    assert tiling.max_diff(map2, map2o) == 0
    print('yes')

if 1:
    print('Cache pointing matrix.', end='\n ...')
    with Timer() as T:
        pix_idx, spin_proj = pe.pointing_matrix(ptg, ofs, None, None)

    print('Forward project using precomputed pointing matrix.',
          end='\n ...')
    pp = so3g.ProjEng_Precomp_Tiled()
    with Timer() as T:
        map1p = pp.to_map(tiling.new_map(3), pix_idx, spin_proj, sig_list, None)

    print('and also weights', end='\n ... ')
    with Timer() as T:
        map2p = pp.to_weight_map(tiling.new_map((3,3)),pix_idx,spin_proj,None,None)

    print('Checking that precomp and on-the-fly forward calcs agree: ',
          end='\n ... ')
    assert tiling.max_diff(map1, map1p) == 0
    assert tiling.max_diff(map2, map2p) == 0
    print('yes')

    print('Reverse projection using precomputed pointing',
          end='\n ...')
    with Timer() as T:
        sig1p = pp.from_map(map1, pix_idx, spin_proj, None, None)

    print('Checking that it agrees with on-the-fly',
          end='\n ...')
    sig1f = pe.from_map(map1, ptg, ofs, None, None)
    thresh = max([m.std()*1e-6 for m in map1[1] if m is not None])
    assert max([np.abs(a - b).max() for a, b in zip(sig1f, sig1p)]) < thresh
    print('yes')

print('Plotting...')
import pylab as pl
gs1 = pl.matplotlib.gridspec.GridSpec(*tiling.tilage)

for y in range(tiling.tilage[0]):
    for x in range(tiling.tilage[1]):
        m = map1[1][y*tiling.tilage[1] + x]
        if m is not None:
            ax = pl.subplot(gs1[y, x])
            ax.imshow(m[0], cmap='gray')

#ax = pl.subplot(gs1[1,:])
#ax.plot(sig1[0])
pl.show()

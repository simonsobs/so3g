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
parser.add_argument('--tiled', action='store_true')
args = parser.parse_args()
system = args.system
print('Using system: %s' % system)

# Map space -- args for Pixelization.
pxz = (300,250,0.00005,0.00005,0.,0.)

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

pe = test_utils.get_proj(system, 'TQU', pxz, tiled=args.tiled)
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
    #pix = np.empty(sig.shape[1:] + (pxz.index_count,), 'int32')
    pix = np.empty(sig.shape[1:] + (pe.index_count,), 'int32')
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
    print('From map into time-domain (TQU)', end='\n ... ')
    map1 = pe.zeros(3) + np.array([1,0,0])[:,None,None]
    with Timer() as T:
        sig1 = pe.from_map(map1, ptg, ofs, None, None)

if 1:
    print('From time into map-domain (TQU)', end='\n ... ')
    if args.tiled:
        map0 = pxz.zeros(3)
        map0 = [(map0.shape, map0.shape), [map0]]
    else:
        map0 = pe.zeros(3)
    sig_list = [x for x in sig[0]]
    with Timer() as T:
        map1 = pe.to_map(map0,ptg,ofs,sig_list,None)

if 1:
    print('Map-domain again but with None for input map', end='\n ... ')
    with Timer() as T:
        map1 = pe.to_map(None,ptg,ofs,sig_list,None)

if 1:
    print('Forward project weights (TQU)', end='\n ... ')
    map0 = pe.zeros((3, 3))
    with Timer() as T:
        map2 = pe.to_weight_map(map0,ptg,ofs,None,None)

if 1:
    print('Weights again but with None for input map', end='\n ... ')
    with Timer() as T:
        map2 = pe.to_weight_map(None,ptg,ofs,None,None)

stop_omp

print('Compute pixel_ranges (OMP prep)... ', end='\n ... ')
with Timer():
    Ivals = pe.pixel_ranges(ptg, ofs)

if 1:
    print('Forward projection (TQU) with OMP (%s): ' % n_omp, end='\n ... ')
    with Timer() as T:
        map1o = pe.to_map_omp(None,ptg,ofs,sig_list,None,Ivals)

if 1:
    print('Reverse projection (TQU)', end='\n ... ')
    sig1 = [x for x in np.zeros((n_det,n_t), 'float32')]
    sig1 = np.random.uniform(size=(n_det,n_t)).astype('float32')
    with Timer() as T:
        #sig1 =
        pe.from_map(map1, ptg, ofs, sig1, None)

if 1:
    print('Forward project weights (TQU) with OMP (%s): ' % n_omp, end='\n ... ')
    with Timer() as T:
        map2o = pe.to_weight_map_omp(None,ptg,ofs,None,None,Ivals)

    print('Checking that OMP and non-OMP forward calcs agree: ', end='\n ... ')
    assert abs(map1 - map1o).max() == 0
    assert abs(map2 - map2o).max() == 0
    print('yes')


if 1:
    print('Cache pointing matrix.', end='\n ...')
    with Timer() as T:
        pix_idx, spin_proj = pe.pointing_matrix(ptg, ofs, None, None)

    print('Forward project using precomputed pointing matrix.',
          end='\n ...')
    pp = so3g.ProjEng_Precomp_()
    map1p = pxz.zeros(3)
    with Timer() as T:
        pp.to_map(map1p, pix_idx, spin_proj, sig_list, None)

    print('and also weights', end='\n ... ')
    map2p = np.zeros((3,3) + pxz.zeros(1).shape[1:])
    with Timer() as T:
        map2p = pp.to_weight_map(map2p,pix_idx,spin_proj,None,None)

    print('Checking that precomp and on-the-fly forward calcs agree: ',
          end='\n ... ')
    assert abs(map1 - map1p).max() == 0
    assert abs(map2 - map2p).max() == 0
    print('yes')

    print('Reverse projection using precomputed pointing',
          end='\n ...')
    with Timer() as T:
        sig1p = pp.from_map(map1, pix_idx, spin_proj, None, None)

    print('Checking that it agrees with on-the-fly',
          end='\n ...')
    sig1f = pe.from_map(map1, ptg, ofs, None, None)
    thresh = map1[0].std() * 1e-6
    assert max([np.abs(a - b).max() for a, b in zip(sig1f, sig1p)]) < thresh
    print('yes')

print('Plotting...')
import pylab as pl
gs1 = pl.matplotlib.gridspec.GridSpec(2, 3)

for axi in range(3):
    ax = pl.subplot(gs1[0,axi])
    ax.imshow(map1[axi], cmap='gray')
    ax.set_title('TQU'[axi])

ax = pl.subplot(gs1[1,:])
ax.plot(sig1[0])
pl.show()

import so3g
import numpy as np
import pylab as pl

from test_utils import Timer, Qmul, Qroti
import os


# Command line option(s)...
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--system', '-s', choices=[
    'flat',
    'qcyl',
    'qzen',
], default='flat')
args = parser.parse_args()
system = args.system
print('Using system: %s' % system)

# Map space.
pxz = so3g.Pixelizor2_Flat(300,250,0.00005,0.00005,0,0)

# Samples
n_det = 2000
n_t = 25000

# Boresight
phi = np.arange(n_t) / n_t * 6.28 * 3
r = np.arange(n_t) / n_t * .006
x = r*np.cos(phi) + .004 + .0005 * np.random.uniform(size=len(phi))
y = r*np.sin(phi) + .008 + .0005 * np.random.uniform(size=len(phi))

# offsets
dr = .002
polphi = 6.28 * np.arange(n_det) / n_det
dx, dy = dr * np.cos(polphi), dr * np.sin(polphi)

if system == 'flat':
    pe = so3g.ProjEng_Flat_TQU(pxz)
    # At each time step, boresight is (x, y, cos(phi), sin(phi))
    ptg = np.zeros((n_t, 4))
    ptg[...,0] = x
    ptg[...,1] = y
    ptg[...,2] = 1.
    ptg[...,3] = 0.
    ofs = np.transpose([dx, dy, np.cos(polphi), np.sin(polphi)])

elif system == 'qcyl':
    pe = so3g.ProjEng_QCyl_TQU(pxz)
    # boresight needs to point to equinox...
    ptg = Qmul(Qroti(2, x),
               Qroti(1, np.pi/2 - y),
               Qroti(2, np.pi))
    ofs = Qmul(Qroti(0, dx),
               Qroti(1,-dy),
               Qroti(2, polphi))

elif system == 'qzen':
    pe = so3g.ProjEng_QZen_TQU(pxz)
    # boresight needs to point to pole...
    ptg = Qmul(Qroti(1, y),
               Qroti(0, x))
    ofs = Qmul(Qroti(0, dx),
               Qroti(1,-dy),
               Qroti(2, polphi))

sig = np.ones((n_det,n_t)) * .5

print('Note the problem size is %i x %i = %.3fM samples.\n' %
      (n_det, n_t, n_det*n_t/1e6))

if 1:
    n_omp = int(os.getenv('OMP_NUM_THREADS'))

    print('Compute pixel_ranges (OMP prep)... ', end='')
    with Timer():
        Ivals = pe.pixel_ranges(ptg, ofs)
        
    # This also works:
    if 0:
      with Timer():
        print('Generating thread-striped map.')
        
        map0 = pxz.zeros(3)
        n_split_ax = map0.shape[-2]
        for i in range(n_omp):
            map0[0,i*n_split_ax//n_omp:(i+1)*n_split_ax//n_omp] = i+1

        print('Project into timestream.')
        stripes = sig * 0 - 1
        pe.from_map(map0, ptg, ofs, stripes, None)

        # Make full set of intervals - n-
        print('Make intervals (%i OMP)' % n_omp)
        Ivals = []
        for i in range(n_omp):
            Ivals.append([
                so3g.IntervalsInt32.from_mask(
                    (st==i).astype('int8'), 1)[0]
                for st in stripes[0]])

    print('Project to map.')
    print('No OMP: ', end='')
    with Timer() as T:
        map1 = pe.to_map(None,ptg,ofs,sig,None)

    print('OMP (%i): ' % n_omp, end='')
    with Timer() as T:
        map1o = pe.to_map_omp(None,ptg,ofs,sig,None,Ivals)

    assert(np.all(map1 == map1o))

    pl.imshow(map1[0])
    pl.show()
    
    print('Reverse projection (TQU)', end='\n ... ')
    with Timer() as T:
        sig1 = pe.from_map(map1, ptg, ofs, None, None)

if 1:
    print('Forward project weights (TQU)')
    map0 = None #pxz.zeros(3)
    print('No OMP: ', end='')
    with Timer() as T:
        map2 = pe.to_weight_map(None,ptg,ofs,None,None)

    print('OMP (%i): ' % n_omp, end='')
    with Timer() as T:
        map2o = pe.to_weight_map_omp(None,ptg,ofs,None,None,Ivals)

    assert(np.all(map2 == map2o))

        
print('Plotting...')
import pylab as pl
gs1 = pl.matplotlib.gridspec.GridSpec(2, 3)

for axi in range(3):
    ax = pl.subplot(gs1[0,axi])
    ax.imshow(map1[axi], cmap='gray')
    ax.set_title('TQU'[axi])

ax = pl.subplot(gs1[1,:])
ax.plot(sig1[500])
pl.show()

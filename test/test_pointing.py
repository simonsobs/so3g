import so3g
import so3g.proj as SP
from spt3g import core
import numpy as np
import pylab as pl

from pixell import enmap

from test_utils import Timer, Qmul, Qroti

# Map space.
pxz = so3g.Pixelizor2_Flat(300,250,0.00005,0.00005,0,0,0,0)

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

# At each time step, boresight is (x, y, cos(phi), sin(phi))
system = 'qcyl'

if system == 'flat':
    pe = so3g.ProjectionEngine2(pxz)
    ptg = np.zeros((n_t, 4))
    ptg[...,0] = x
    ptg[...,1] = y
    ptg[...,2] = 1.
    ptg[...,3] = 0.
    ofs = np.transpose([dx, dy, np.cos(polphi), np.sin(polphi)])

elif system == 'qcyl':
    pe = so3g.ProjectionEngine2QC(pxz)
    # boresight needs to point to equinox...
    ptg = Qmul(Qroti(2, x),
               Qroti(1, np.pi/2 - y),
               Qroti(2, np.pi))
    ofs = Qmul(Qroti(0, dx),
               Qroti(1,-dy),
               Qroti(2, polphi))

elif system == 'qzen':
    pe = so3g.ProjectionEngine2QZ(pxz)
    # boresight needs to point to pole...
    ptg = Qmul(Qroti(1, y),
               Qroti(0, x))
    ofs = Qmul(Qroti(0, dx),
               Qroti(1,-dy),
               Qroti(2, polphi))

sig = np.ones((1,n_det,n_t)) * .5

map0 = pxz.zeros(1)

if 1:
    coo = np.empty(sig.shape[1:] + (4,), 'double')
    print('Compute and return coordinates only.', end='\n ... ')
    with Timer() as T:
        pe.coords(ptg,ofs[:,:],coo)

    pl.plot(coo[0,:,0],
            coo[0,:,1])
    pl.show()

    del coo

    print('Compute coords and pixels and return pixels.', end='\n ... ')
    pix = np.empty(sig.shape[1:], 'int32')
    with Timer() as T:
        pe.pixels(ptg,ofs,pix)

    del pix

if 1:
    print('Forward projection (TQU)', end='\n ... ')
    map0 = None #pxz.zeros(3)
    with Timer() as T:
        map1 = pe.to_map(None,ptg,ofs,sig,None)

    print('Reverse projection (TQU)', end='\n ... ')
    with Timer() as T:
        sig1 = pe.from_map(map1, ptg, ofs, None, None)

if 1:
    print('Forward project weights (TQU)', end='\n ... ')
    map0 = None #pxz.zeros(3)
    with Timer() as T:
        map2 = pe.to_weight_map(None,ptg,ofs,None,None)

print('Plotting...')
import pylab as pl
gs1 = pl.matplotlib.gridspec.GridSpec(2, 3)

for axi in range(3):
    ax = pl.subplot(gs1[0,axi])
    ax.imshow(map1[axi], cmap='gray')
    ax.set_title('TQU'[axi])

ax = pl.subplot(gs1[1,:])
ax.plot(sig1[0,500])
pl.show()

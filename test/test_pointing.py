import so3g
import so3g.proj as SP
from spt3g import core
import numpy as np

from pixell import enmap

from test_utils import Timer

pxz = so3g.Pixelizor(300,250,0.00005,0.00005,0,0,0,0)
pe = so3g.ProjectionEngine0(pxz)
map0 = pxz.zeros(1)

n_det = 2000
n_t = 25000
ptg = np.zeros((n_t, 4))
phi = np.arange(n_t) / n_t * 6.28 * 3
r = np.arange(n_t) / n_t * .006
ptg[...,0] = r*np.cos(phi) + .004 + .0005 * np.random.uniform(size=len(phi))
ptg[...,1] = r*np.sin(phi) + .008 + .0005 * np.random.uniform(size=len(phi))
ptg[...,2] = 1.
ptg[...,3] = 0.

ophi = 6.28 * np.arange(n_det) / n_det
ofs = np.transpose([np.cos(ophi), np.sin(ophi), np.cos(ophi), np.sin(ophi)])
ofs[:,:2] *= 0.001

sig = np.ones((1,n_det,n_t))

#with Timer() as T:
#    map1 = pe.to_map(map0,ptg,ofs,sig,wts)
#
#print('%i out of %i' % (map1.sum(), n_det*n_t))
#
pe = so3g.ProjectionEngine2(pxz)

coo = np.empty(sig.shape[1:] + (4,), 'double')
print('Compute and return coordinates only.', end='\n ... ')
with Timer() as T:
    pe.coords(ptg,ofs[:,:],coo)

del coo

print('Compute coords and pixels and return pixels.', end='\n ... ')
pix = np.empty(sig.shape[1:], 'int32')
with Timer() as T:
    pe.pixels(ptg,ofs,pix)

del pix

print('Forward projection (TQU)', end='\n ... ')
with Timer() as T:
    map1 = pe.to_map(None,ptg,ofs,sig,None)

print('Reverse projection (TQU)', end='\n ... ')
sig[:] = 0
with Timer() as T:
    pe.from_map(map1, ptg, ofs, sig, None)

print('Plotting...')
import pylab as pl
gs1 = pl.matplotlib.gridspec.GridSpec(2, 3)

for axi in range(3):
    ax = pl.subplot(gs1[0,axi])
    ax.imshow(map1[axi], cmap='gray')
    ax.set_title('TQU'[axi])

ax = pl.subplot(gs1[1,:])
ax.plot(sig[0,0])
pl.show()

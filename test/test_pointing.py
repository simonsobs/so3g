import so3g
import so3g.proj as SP
from spt3g import core
import numpy as np
import pylab as pl

from pixell import enmap

from test_utils import Timer

pe = so3g.ProjectionEngine()
map0 = pe.zeros(None)[:,:,None]

n_det = 2000
n_t = 25000
ptg = np.zeros((n_t, 4))
phi = np.arange(n_t) / n_t * 6.28 * 3
r = np.arange(n_t) / n_t * .006
ptg[...,0] = r*np.cos(phi) + .004 + .0005 * np.random.uniform(size=len(phi))
ptg[...,1] = r*np.sin(phi) + .008 + .0005 * np.random.uniform(size=len(phi))
ptg[...,2] = np.cos(phi)
ptg[...,3] = np.sin(phi)

ophi = 6.28 * np.arange(n_det) / n_det
ofs = .001 * np.transpose([np.cos(ophi), np.sin(ophi)])

wts = None
sig = np.ones((1,n_det,n_t))

with Timer() as T:
    map1 = pe.to_map(map0,ptg,ofs,sig,wts)

print('%i out of %i' % (map1.sum(), n_det*n_t))

pe = so3g.ProjectionEngine2()
map0 = pe.zeros(None)
map0 = np.zeros(map0.shape + (3,), map0.dtype)
#map0 = np.zeros((3,) + map0.shape, map0.dtype).transpose((1,2,0))
print(map0.shape)

with Timer() as T:
    map1 = pe.to_map(map0,ptg,ofs,sig,wts)

for axi,ax in enumerate(pl.subplots(1,3)[1]):
    ax.imshow(map1[...,axi])
    ax.set_title('TQU'[axi])

pl.show()

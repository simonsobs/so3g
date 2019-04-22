import so3g
import so3g.proj as SP
from spt3g import core
import numpy as np
import pylab as pl

from pixell import enmap

from test_utils import Timer

x = so3g.GnomonicGridder()
map0 = x.zeros(None)[:,:,None]

#n_det = 200
#n_t = 2500
#ptg = np.zeros((n_det, n_t, 2))
#phi = np.arange(n_t) / n_t * 6.28 * 3
#r = np.arange(n_t) / n_t * .006
#ptg[...,0] = r*np.cos(phi) + .004 + .0005 * np.random.uniform(size=len(phi))
#ptg[...,1] = r*np.sin(phi) + .008 + .0005 * np.random.uniform(size=len(phi))
#
#
#wts = np.zeros((1,n_det,1))
#wts[:] = 1.
#sig = np.ones((1,n_det,n_t))
#
#map1 = x.to_map(map0,ptg,sig,wts)
#
#pl.imshow(map1[...,0])
#pl.show()
#

n_det = 2000
n_t = 25000
ptg = np.zeros((n_t, 2))
phi = np.arange(n_t) / n_t * 6.28 * 3
r = np.arange(n_t) / n_t * .006
ptg[...,0] = r*np.cos(phi) + .004 + .000 * np.random.uniform(size=len(phi))
ptg[...,1] = r*np.sin(phi) + .008 + .000 * np.random.uniform(size=len(phi))

ophi = 6.28 * np.arange(n_det) / n_det
ofs = .001 * np.transpose([np.cos(ophi), np.sin(ophi)])

#wts = np.zeros((1,n_det,1))
#wts[:] = 0.5
wts = None
sig = np.ones((1,n_det,n_t))

with Timer() as T:
    map1 = x.to_map2(map0,ptg,ofs,sig,wts)

print(map1.sum())
pl.imshow(map1[...,0])
pl.show()

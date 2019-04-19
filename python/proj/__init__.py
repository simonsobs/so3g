import so3g
from spt3g import core

import numpy as np

"""We're typically converting from some set of vectors of time-ordered
information (e.g. boresight az-el) to per-detector pointing
(ra,dec,gamma) for each detector.

Perturbation approaches allow us to compute the boresight pointing,
and some derivatives, and then use detector offsets.

First order corrections to detector offsets (a.k.a. 2nd order
perturbations) arise due to non-linear effects such as:

- atmospheric refraction, which depends on elevation
- aberration, which depends on (in some sense) RA, dec
"""
import time

def test():
    # Create a map.
    from pixell import enmap
    map0 = enmap.zeros((1024,1024))
    f = core.G3Frame()
    f['map_out'] = map0
    
    

def profile():
    import healpy
    N = 1000 * 400 * 60 * 10
    print('Create vector (%i).' % N)
    z0 = np.random.uniform(size=N)
    times = [time.time()]
    for op in ['add', 'scale', lambda x: healpy.ang2pix(2048, x, .3), np.cos, np.arcsin, np.sqrt, ]:
        print(op)
        if str(op) == 'add':
            a = z0 + 10.
        elif str(op) == 'scale':
            a = z0 * 10.
        else:
            a = op(z0)
        times.append(time.time())
        print(times[-1] - times[-2])
    

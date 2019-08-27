import time
import numpy as np
import collections

import so3g

class Timer(): 
    def __init__(self, block_name=None): 
        self.t0 = time.time()
          
    def __enter__(self): 
        return self
      
    def report(self):
        return time.time() - self.t0

    def __exit__(self, exc_type, exc_value, exc_traceback): 
        print('Timer block exits after %.6f seconds' % (time.time() - self.t0))


def Qmul(q1, q2, *args):
    qout = q1[...,:1] * q2
    qout[...,1:] += q1[...,1:] * q2[...,:1]
    qout[...,0] -= (q1[...,1:] * q2[...,1:]).sum(axis=-1)
    qout[...,1] += q1[...,2]*q2[...,3] - q1[...,3]*q2[...,2]
    qout[...,2] += q1[...,3]*q2[...,1] - q1[...,1]*q2[...,3]
    qout[...,3] += q1[...,1]*q2[...,2] - q1[...,2]*q2[...,1]
    if len(args) == 0:
        return qout
    return Qmul(qout, args[0], *args[1:])

def Qroti(n,phi):
    """
    Euler quaternion, appropriate for rotating by angle phi about axis
    n=(0,1,2).
    """
    phi = np.asarray(phi)/2
    out = np.zeros(np.asarray(phi).shape + (4,))
    out[...,0  ] = np.cos(phi)
    out[...,n+1] = np.sin(phi)
    return out

proj_dict = collections.OrderedDict([
    ('flat', 'Flat'),
    ('car' , 'CAR'),
    ('cea' , 'CEA'),
    ('arc' , 'ARC'),
    ('tan' , 'TAN'),
    ('zea' , 'ZEA'),
    ])

def get_proj(coord_sys, pol_sys, pxz=None):
    assert pol_sys in ['T', 'TQU', 'QU']
    name = 'ProjEng_{}_{}'.format(proj_dict[coord_sys], pol_sys)
    cls = getattr(so3g, name) # ProjEng_X_Y
    if pxz is None:
        return cls
    return cls(pxz)

def get_boresight_quat(system, x, y, gamma=None):
    if gamma is None:
        gamma = 0

    if system == 'flat':
        # At each time step, boresight is (x, y, cos(phi), sin(phi))
        n_t = len(x)
        ptg = np.zeros((n_t, 4))
        ptg[...,0] = x
        ptg[...,1] = y
        ptg[...,2] = np.cos(gamma)
        ptg[...,3] = np.sin(gamma)

    elif system in ['car', 'cea']:
        # boresight needs to point to equinox...
        ptg = Qmul(Qroti(2, x),
                   Qroti(1, np.pi/2 - y),
                   Qroti(2, np.pi - gamma))

    elif system in['arc', 'tan', 'zea']:
        # boresight needs to point to pole...
        ptg = Qmul(Qroti(1, y),
                   Qroti(0, x),
                   Qroti(2, gamma))

    else:
        raise ValueError('Unknown system: "%s"' % system)

    return ptg

def get_offsets_quat(system, dx, dy, polphi):
    if system == 'flat':
        ofs = np.transpose([dx, dy, np.cos(polphi), np.sin(polphi)])
    else:
        ofs = Qmul(Qroti(0, dx),
                   Qroti(1,-dy),
                   Qroti(2, polphi))
    return ofs

def linalg_pinv(wmap):
    try:
        iwmap = np.linalg.pinv(wmap)
    except ValueError:
        # np.linalg.pinv only runs on stacks for numpy>=1.14 (Jan 2018),
        # so we hack an alternative.
        iwmap = np.empty(wmap.shape, wmap.dtype)
        for i in range(wmap.shape[0]):
            for j in range(wmap.shape[1]):
                iwmap[i,j] = np.linalg.pinv(wmap[i,j])
    return iwmap

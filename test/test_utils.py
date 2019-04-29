import time
import numpy as np

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

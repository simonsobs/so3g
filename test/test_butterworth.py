import so3g
import numpy as np

import ctypes

bank = so3g.BFilterBank() \
           .add(so3g.BFilterParams(32092,15750,14,3,5)) \
           .add(so3g.BFilterParams(31238,14895,14,3,12)) \
           .init(1)

n = 100
a1 = np.ones(n, 'int32')
b1 = 0*a1

bank.apply(a1, b1)
assert np.all(b1[::10] == [   0,    5,   45,  151,  329,  560,  809, 1042, 1230, 1357])

a2 = a1.astype('float32')
b2 = a2*0

bank.init(1)
bank.apply(a2, b2)

assert np.all(b1 == b2)

import so3g
core = so3g.spt3g_core

import numpy as np

import gc
print('Creating y')
y = np.arange(100,200000, dtype=float)
x = so3g.G3Ndarray(y)
print(x)

# write that out?

w = core.G3Writer('out.g3')
f = core.G3Frame()
f['x'] = x
w.Process(f)
del w

# read that back?

print('Readback:')
for f in core.G3File('out.g3'):
    print(f['x'])


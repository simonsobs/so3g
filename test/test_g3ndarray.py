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

# And the autoconverter?
f = so3g.SOFrame()
f['y'] = 1
f['x'] = np.arange(1000.)
w.Process(f)

del w

# read those back?

print('Readback:')
for f in core.G3File('out.g3'):
    print(f['x'])


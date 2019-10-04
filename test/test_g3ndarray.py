import so3g
core = so3g.spt3g_core

import numpy as np
import gc
import unittest

# Don't require pixell for testing
try:
    from pixell import enmap
    pixell_found = True
except ModuleNotFoundError:
    pixell_found = False

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
f = core.G3Frame()
f['y'] = 1
f['x'] = np.arange(1000.)
w.Process(f)

del w

# read those back?

print('Readback:')
for f in core.G3File('out.g3'):
    print(f['x'])

class TestG3NDArray(unittest.TestCase):
    """TestCase for testing G3NDArrays."""
    @unittest.skipIf(pixell_found is False, "pixell not found")
    def test_enmap(self):
        """Test writing enmap objects to G3Frames."""
        # enmap?
        m = enmap.zeros((100,20))
        f = core.G3Frame()
        f['test_map'] = m

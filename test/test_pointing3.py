import so3g
import numpy as np
import pylab as pl

# Test pointing support functions.
#### q test.
from so3g.proj import quat

iso_angles    = (.1, .2, .05)
lonlat_angles = (iso_angles[1], np.pi/2 - iso_angles[0], iso_angles[2])

# Create rotation both ways.
q0 = quat.rotation_iso(*iso_angles)
q1 = quat.rotation_lonlat(*lonlat_angles)

# Decompose to recover originals.
iso_angles1    = quat.decompose_iso(q0)
lonlat_angles1 = quat.decompose_lonlat(q0)

print('These quats should be the same:\n {}\n {}\n'.format(q0, q1))

print('These angle sets should be very close:\n {}\n {}\n'\
      .format(iso_angles, iso_angles1))

print('These angle sets should be very close:\n {}\n {}\n'\
      .format(lonlat_angles, lonlat_angles1))


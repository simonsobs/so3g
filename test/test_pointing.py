import so3g
import so3g.proj as SP
from spt3g import core
import numpy as np

from pixell import enmap


x = so3g.GnomonicGridder()
map0 = x.zeros(None)
ptg = np.zeros((2,100))
ptg[0,:] = .00033
ptg[1,:] = .00044

map1 = x.to_map(map0,ptg,None,None)


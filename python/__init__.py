# For our compiled libraries to load, the spt3g.core library must already be loaded.
from spt3g import core as spt3g_core

# But we use our own load_pybindings.  And our library is called
# libso3g.so, not so3g.so.  How to I feel about that?  libso-so.
from .load_pybindings import load_pybindings
load_pybindings([__path__[0] + '/libso3g'], name='so3g')

__version__ = version()

# And python stuff.
from . import hkagg
from . import coords

import numpy as np
from .soframe import SOFrame
SOFrame.setitem_converters[np.ndarray]  = lambda a: G3Ndarray(a)
SOFrame.getitem_converters[G3Ndarray] = lambda a: a.to_array()

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

from .soframe import SOFrame
# Numpy arrays in frames
import numpy as np
SOFrame.setitem_converters[np.ndarray] = lambda a: G3Ndarray(a)
SOFrame.getitem_converters[G3Ndarray]  = lambda a: a.to_array()
# Astropy wcs in frames
import astropy.wcs
SOFrame.setitem_converters[astropy.wcs.WCS] = lambda a: G3WCS(a.to_header_string())
SOFrame.getitem_converters[G3WCS]           = lambda a: astropy.wcs.WCS(a.header)
# Enmaps
from pixell import enmap
SOFrame.setitem_converters[enmap.ndmap] = lambda a: G3Ndmap(a, a.wcs.to_header_string())
SOFrame.getitem_converters[G3Ndmap]     = lambda a: enmap.ndmap(a.data.to_array(), astropy.wcs.WCS(a.wcs.header))

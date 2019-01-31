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

from spt3g.core import G3Frame
G3Frame.__doc__ = "Monkey patched in so3g.__init__ to add support for numpy arrays etc.\n\n" + G3Frame.__doc__
orig_getitem = G3Frame.__getitem__
orig_setitem = G3Frame.__setitem__

G3Frame.getitem_converters = {}
G3Frame.setitem_converters = {}

def patched_getitem(self, key):
	val = orig_getitem(self, key)
	try: val = G3Frame.getitem_converters[type(val)](val)
	except KeyError: pass
	return val
def patched_setitem(self, key, val):
	try: val = G3Frame.setitem_converters[type(val)](val)
	except KeyError: pass
	orig_setitem(self, key, val)

G3Frame.__getitem__ = patched_getitem
G3Frame.__setitem__ = patched_setitem

# Numpy arrays in frames
import numpy as np
G3Frame.setitem_converters[np.ndarray] = lambda a: G3Ndarray(a)
G3Frame.getitem_converters[G3Ndarray]  = lambda a: a.to_array()
# Astropy wcs in frames
import astropy.wcs
G3Frame.setitem_converters[astropy.wcs.WCS] = lambda a: G3WCS(a.to_header_string())
G3Frame.getitem_converters[G3WCS]           = lambda a: astropy.wcs.WCS(a.header)
# Enmaps
from pixell import enmap
G3Frame.setitem_converters[enmap.ndmap] = lambda a: G3Ndmap(a, a.wcs.to_header_string())
G3Frame.getitem_converters[G3Ndmap]     = lambda a: enmap.ndmap(a.data.to_array(), astropy.wcs.WCS(a.wcs.header))

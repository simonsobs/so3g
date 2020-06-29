"""The code in this file helps to extend the automatic conversions
between G3 types and python types when getting/setting objects a
G3Frame.  This permits seamless assignment and retrieval of, for
example, numpy arrays, which will automatically be encoded to or
decoded from the C++ type, G3Ndarray.

For example:

      frame = core.G3Frame()
      frame["data"] = np.arange(10.)
      print(frame["data"])

This is "monkey patching" because it changes the behavior of objects
originally defined in spt3g library.  The behavior can be disabled
with the config variable patch_g3frame.
"""

import so3g
from spt3g.core import G3Frame

orig_getitem = None
orig_setitem = None


def patched_getitem(self, key):
    val = orig_getitem(self, key)
    try:
        val = G3Frame.getitem_converters[type(val)](val)
    except KeyError:
        pass
    return val


def patched_setitem(self, key, val):
    try:
        val = G3Frame.setitem_converters[type(val)](val)
    except KeyError:
        pass
    orig_setitem(self, key, val)


def _try_import(module_name, insistance):
    """Maybe import a module and maybe raise an error if the import fails.
    """
    if insistance in [False, 'no', 'false', 'False']:
        return None
    module = None
    try:
        module = __import__(module_name)
    except ImportError as e:
        if insistance in [True, 'yes', 'true', 'True']:
            raise(e)
    return module


def set_frame_hooks(config={}):
    if not config.get('patch_g3frame', True):
        return

    G3Frame.__doc__ = ("Monkey patched by so3g to add support for numpy "
                       "arrays etc.\n\n" + G3Frame.__doc__)

    # Do not trash orig_* if we've been here before.
    if not hasattr(G3Frame, 'getitem_converters'):
        global orig_getitem, orig_setitem
        orig_getitem = G3Frame.__getitem__
        orig_setitem = G3Frame.__setitem__

    # But do reset the converters.
    G3Frame.getitem_converters = {}
    G3Frame.setitem_converters = {}

    G3Frame.__getitem__ = patched_getitem
    G3Frame.__setitem__ = patched_setitem

    # Always do numpy.
    import numpy as np
    # Numpy arrays in frames
    G3Frame.setitem_converters[np.ndarray] = lambda a: so3g.G3Ndarray(a)
    G3Frame.getitem_converters[so3g.G3Ndarray] = lambda a: a.to_array()

    has_astropy = False
    use_astropy = config.get('use_astropy', 'try')
    astropy = _try_import('astropy.wcs', use_astropy)
    if astropy is not None:
        has_astropy = True
        G3Frame.setitem_converters[astropy.wcs.WCS] = \
            lambda a: so3g.G3WCS(a.to_header_string())
        G3Frame.getitem_converters[so3g.G3WCS] = \
            lambda a: astropy.wcs.WCS(a.header)

    use_pixell = config.get('use_pixell', 'try')
    pixell = _try_import('pixell.enmap', use_pixell)
    if pixell is not None and has_astropy:
        G3Frame.setitem_converters[pixell.enmap.ndmap] = \
            lambda a: so3g.G3Ndmap(a, a.wcs.to_header_string())
        G3Frame.getitem_converters[so3g.G3Ndmap] = \
            lambda a: pixell.enmap.ndmap(a.data.to_array(),
                                         astropy.wcs.WCS(a.wcs.header))

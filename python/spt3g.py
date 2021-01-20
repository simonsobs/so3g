"""
This package simply provides a universal way to import spt3g, either from a
bundled subpackage or from somewhere on the filesystem.
"""

try:
    from .spt3g_internal import core, dfmux, gcp, calibration, maps
except:
    # Not bundled
    try:
        from spt3g import core, dfmux, gcp, calibration, maps
    except:
        raise ImportError("Cannot import either the internal or external spt3g!")

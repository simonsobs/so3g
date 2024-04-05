"""
This package simply provides a universal way to import spt3g, either from a
bundled subpackage or from somewhere on the filesystem.
"""

import sys

try:
    from . import spt3g_internal
    sys.modules["spt3g"] = sys.modules["so3g.spt3g_internal"]
    sys.modules["spt3g"].__name__ = "spt3g"
    del sys.modules["so3g.spt3g_internal"]
    from spt3g import core, __version__, __file__
except:
    # Not bundled
    try:
        from spt3g import core, __version__, __file__
    except:
        raise ImportError("Cannot import either the internal or external spt3g!")

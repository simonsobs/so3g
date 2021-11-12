import os
import sys

if os.getenv('DOCS_BUILD') == '1':
    from ._libso3g_docstring_shells import *
else:
    # For our compiled libraries to load, the spt3g.core library must already be loaded.
    from . import spt3g
    from .spt3g import core as spt3g_core

    # Forcibly export our spt3g wrapper module as the global spt3g import
    sys.modules["spt3g"] = spt3g

    # FIXME:  We could restrict this list to just the "top-level" symbols in the
    # the compiled extension, but that might break external packages that use, for
    # example, "so3g.HKFrameType" instead of "so3g.hk.HKFrameType".  For now, symbols
    # that are specific to the hk and proj submodules are available in both places.
    from .libso3g import *

# Version is computed by versioneer.
__version__ = version()

if os.getenv('DOCS_BUILD') != '1':
    # Instance configuration.
    from .config import get_config
    instance_config = get_config()
    del get_config

    # (Possibly) monkey patch the G3Frame object with hooks for
    # so3g data types.
    from .soframe import set_frame_hooks
    set_frame_hooks(instance_config)
    del set_frame_hooks

# Other python modules.
from . import hk
from . import proj

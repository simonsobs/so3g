import os

if os.getenv('DOCS_BUILD') == '1':
    from ._libso3g_docstring_shells import *
else:
    # For our compiled libraries to load, the spt3g.core library must already be loaded.
    from . import spt3g
    from spt3g import core as spt3g_core

    # Our library is called libso3g.{suffix}, but will load into module
    # namespace so3g.
    from .load_pybindings import load_pybindings
    load_pybindings([__path__[0] + '/libso3g'], name='so3g')

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

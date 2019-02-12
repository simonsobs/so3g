# For our compiled libraries to load, the spt3g.core library must already be loaded.
from spt3g import core as spt3g_core

# Our library is called libso3g.{suffix}, but will load into module
# namespace so3g.
from .load_pybindings import load_pybindings
load_pybindings([__path__[0] + '/libso3g'], name='so3g')

# Version is computed by versioneer.
__version__ = version()

# Instance configuration.
from .config import get_config
instance_config = get_config()
del get_config

# Get the SOFrame object, and possibly enable transparent get/set of
# our G3-compatible objects.
from .soframe import SOFrame, set_frame_hooks
set_frame_hooks(instance_config)
del set_frame_hooks

# Other python modules.
from . import hkagg
from . import coords


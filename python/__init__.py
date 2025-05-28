import os
import numpy as np


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

# Other python modules.
from . import hk
from . import proj

from .g3reader_shim import G3IndexedReader

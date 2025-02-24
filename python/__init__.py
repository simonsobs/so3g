import os


if os.getenv('DOCS_BUILD') == '1':
    from ._libso3g_docstring_shells import *
else:
    # For our compiled libraries to load, the spt3g.core library must already be loaded.
    from . import spt3g
    from spt3g import core as spt3g_core

    # Load all symbols from our compiled extension.  This is needed for backwards
    # compatibility with external packages that expect those symbols to exist based
    # on how the previous loader worked.
    from ._libso3g import *

# Version is computed by versioneer.
__version__ = version()

# Other python modules.
from . import hk
from . import proj

from .g3reader_shim import G3IndexedReader

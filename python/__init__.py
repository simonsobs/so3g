import os
import numpy as np

from .libso3g import version

# Version is defined in the compiled extension.
__version__ = version()

# Other python modules.
from . import hk
#from . import proj

#from .g3reader_shim import G3IndexedReader

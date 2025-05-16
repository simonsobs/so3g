import os
import numpy as np

# Many downstream packages expect all of the libso3g symbols to be exported
# into the top namespace.
from .libso3g import *

# Version is defined in the compiled extension.
__version__ = version()

# Other python modules.
from . import hk
from . import proj

from .g3reader_shim import G3IndexedReader

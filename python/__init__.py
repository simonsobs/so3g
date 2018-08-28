## This is the spt3g way... and works.
from spt3g.core import load_pybindings
load_pybindings(__name__, [__path__[0] + '/so3g'])
## This is another way, which fails in the presence of G3 stuff.
#from .so3g import *

# And python stuff.
from .w import *

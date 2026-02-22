
import numpy as np

from spt3g import core

from . import quat
from . import util
from . import mapthreads

from .wcs import Projectionist, ProjectionistHealpix, Ranges, RangesMatrix
from .coords import CelestialSightLine, EarthlySite, Assembly, FocalPlane, SITES
from .weather import Weather, weather_factory
from .ranges import Ranges, RangesMatrix


DEG = np.pi/180.

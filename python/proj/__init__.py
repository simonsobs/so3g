import so3g
from spt3g import core

from . import quat
from . import util

from .wcs import Projectionist, ProjectionOmpData
from .coords import CelestialSightLine, EarthlySite, Assembly, FocalPlane
from .weather import Weather, weather_factory

import numpy as np

DEG = np.pi/180.

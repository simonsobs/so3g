import so3g
from . import quat

from collections import OrderedDict

import numpy as np

DEG = np.pi / 180.

class EarthlySite:
    def __init__(self, lon, lat, elev):
        """Arguments in degrees.  Positive longitude is East on the Earth."""
        self.lon, self.lat, self.elev = lon, lat, elev

    @classmethod
    def get_named(cls, name):
        return SITES[name]


SITES = {
    'act': EarthlySite(-67.7876 * DEG,-22.9585 * DEG, 5188.),
}
DEFAULT_SITE = 'act'

ERA_EPOCH = 946684800 + 3600 * 12  # Noon on Jan 1 2000.
ERA_POLY = np.array([6.300387486754831,
                     4.894961212823756]) # Operates on "days since ERA_EPOCH".

class CelestialSightLine:
    """Carries a vector of celestial pointing data.

    Instantiation should result in the pointing quaternion being
    stored in self.Q.

    """
    @staticmethod
    def decode_site(site=None):
        """Convert site argument to a Site object."""
        if site is None:
            site = DEFAULT_SITE
        if isinstance(site, EarthlySite):
            return site
        if site in SITES:
            return SITES[site]
        raise ValueError("Could not decode %s as a Site." % site)

    @classmethod
    def naive_az_el(cls, az, el, t, site=None):
        """Construct a sight line from horizon coordinates, az and el (in
        radians) and time t (unix timestamp).

        This will be off by several arcminutes... but less than a
        degree.

        """
        site = cls.decode_site(site)
        assert isinstance(site, EarthlySite)

        self = cls()

        J = (t - ERA_EPOCH) / 86400
        era = np.polyval(ERA_POLY, J)
        lst = era + site.lon

        self.Q = (
            quat.euler(2, lst) * quat.euler(1, np.pi/2 - site.lat) *
            quat.euler(2, np.pi) *
            quat.euler(2, -az)  * quat.euler(1, np.pi/2 - el) *
            quat.euler(2, np.pi))
        return self

    def coords(self, det_offsets=None, output=None):
        """Get the celestial coordinates of each detector at each time.

        Arguments:
          det_offset: A dict or list or array of detector offset
            tuples.  If each tuple has 2 elements or 3 elements, the
            typles are interpreted as X,Y[,phi] coordinates in the
            conventional way.  If 4 elements, it's interpreted as a
            quaternion.  If this argument is None, then the boresight
            pointing is returned.
          output: An optional structure for holding the results.  For
            that to work, each element of output must support the
            buffer protocol.

        Returns:
          If the det_offset was passed in as a dict, then a dict with the same
          keys is returned.  Otherwise a list is returned.  For each
          detector, the corresponding returned object is an array with
          shape (n_samp, 4) where the four components correspond to
          longitude, latitude, cos(gamma), sin(gamma).

        """
        # Get a projector, in CAR.
        px = so3g.Pixelizor2_Flat(1,1,1.,1.,1.,1.)
        p = so3g.ProjEng_CAR_TQU(px)
        # Pre-process the offsets
        collapse = (det_offsets is None)
        if collapse:
            det_offsets = np.array([[1., 0,0,0]])
            if output is not None:
                output = output[None]
        redict = isinstance(det_offsets, dict)
        if redict:
            keys, det_offsets = zip(*det_offsets.items())
            det_offsets = np.array(det_offsets)
            if isinstance(output, dict):
                output = [output[k] for k in keys]
        if np.shape(det_offsets)[1] == 2:
            # XY
            x, y = det_offsets[:,0], det_offsets[:,1]
            theta = (x**2 + y**2)**0.5
            v = np.array([-y, x]) / theta
            v[np.isnan(v)] = 0.
            det_offsets = np.transpose([np.cos(theta/2),
                                        np.sin(theta/2) * v[0],
                                        np.sin(theta/2) * v[1],
                                        np.zeros(len(theta))])
            det_offsets = quat.cu3g.G3VectorQuat(det_offsets.copy())
        output = p.coords(self.Q, det_offsets, output)
        if redict:
            output = OrderedDict(zip(keys, output))
        if collapse:
            output = output[0]
        return output


class Assembly:
    """This class groups together boresight and detector offset
    information for use by the Projectionist.

    The basic implementation simply attachs some number of detectors
    rigidly to a boresight.  But this abstraction layer could
    facilitate more complex arrangements, eventually.

    """
    def __init__(self, keyed=False, collapse=False):
        self.keyed = keyed
        self.collapse = collapse

    @classmethod
    def attach(cls, sight_line, det_offsets):
        keyed = isinstance(det_offsets, dict)
        self = cls(keyed=keyed)
        self.Q = sight_line.Q
        if self.keyed:
            self.keys = list(det_offsets.keys())
            self.dets = [det_offsets[k] for k in self.keys]
        else:
            self.dets = det_offsets
        return self

    @classmethod
    def for_boresight(cls, sight_line):
        self = cls(collapse=True)
        self.Q = sight_line.Q
        self.dets = [np.array([1.,0,0,0])]
        return self


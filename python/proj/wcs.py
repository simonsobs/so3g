import so3g
from . import quat

import numpy as np

class ProjectionInfo:
    @classmethod
    def for_map(cls, emap):
        self = cls()
        wcs = emap.wcs

        # Figure out the relevant shape in the map.
        ax1, ax2 = wcs.wcs.lng, wcs.wcs.lat
        self.naxis = np.array([emap.shape[ax1], emap.shape[ax2]], dtype=int)

        # Get just the celestial part.
        wcs = wcs.celestial
        proj = [c[-3:] for c in wcs.wcs.ctype]
        assert(proj[0] == proj[1])
        proj_name = proj[0] # Projection name
        coord_names = [c[:-3].rstrip('-') for c in wcs.wcs.ctype]

        self.proj_name = proj_name

        # 1. Determine the native -> celestial rotation.
        #
        # As per FITS-II, the native point (phi0,theta0) is mapped to the
        # celestial point (crval1,crval2).  The 3rd Euler angle is given
        # in LONPOLEa, and is the "native longitude of the celestial
        # pole".
        #
        # This can be done generally, without checking the projection.
        # But for now let's assert and rely on the simple special cases.

        assert(proj_name == 'CAR')

        assert(wcs.wcs.phi0 == 0. and wcs.wcs.theta0 == 0.)
        assert(wcs.wcs.crval[1] == 0.)
        self.q_celestial_to_native = quat.euler(2, -wcs.wcs.crval[0] * quat.DEG)

        # 2. Get the needed pixelization information.
        self.cdelt = np.array(wcs.wcs.cdelt)
        self.crpix = np.array(wcs.wcs.crpix)

        return self

    def get_pixelizor(self):
        # All these casts are required because boost-python doesn't
        # like numpy scalars.
        return so3g.Pixelizor2_Flat(int(self.naxis[1]), int(self.naxis[0]),
                                    float(self.cdelt[1]), float(self.cdelt[0]),
                                    float(self.crpix[1]), float(self.crpix[0]))


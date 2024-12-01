import so3g
from . import quat
from .weather import weather_factory

from collections import OrderedDict

import numpy as np

DEG = np.pi / 180.


class EarthlySite:
    def __init__(self, lon, lat, elev, typical_weather=None):
        """Arguments in degrees E, degrees N, meters above sea level.  The
        typical_weather argument should be an so3g Weather object, or
        None.

        """
        self.lon, self.lat, self.elev = lon, lat, elev
        self.typical_weather = typical_weather

    def ephem_observer(self):
        """Return an ephem.Observer corresponding to this location."""
        import ephem
        site = ephem.Observer()
        site.lat =  str(self.lat)
        site.long = str(self.lon)
        site.elev = self.elev
        return site

    def skyfield_site(self, spice_kernel):
        """Return a skyfield VectorSum for this location on 'earth' (starting
        from wgs84).  The spice_kernel is the thing one might get from
        jpllib.SpiceKernel(emphemeris_filename).

        """
        from skyfield.api import wgs84
        earth = spice_kernel['earth']
        return earth + wgs84.latlon(self.lat, self.lon, self.elev)


def _debabyl(deg, arcmin, arcsec):
    return deg + arcmin/60 + arcsec/3600

SITES = {
    'act': EarthlySite(-67.7876, -22.9585, 5188.,
                       typical_weather=weather_factory('toco')),
    # SO coords taken from SO-SITE-HEF-003A on 2020-06-02; altitude
    # set to same as ACT.
    'so_lat': EarthlySite(-_debabyl(67,47,15.68), -_debabyl(22,57,39.47), 5188.,
                          typical_weather=weather_factory('toco')),
    'so_sat1': EarthlySite(-_debabyl(67,47,18.11), -_debabyl(22,57,36.38), 5188.,
                           typical_weather=weather_factory('toco')),
    'so_sat2': EarthlySite(-_debabyl(67,47,17.28), -_debabyl(22,57,36.35), 5188.,
                           typical_weather=weather_factory('toco')),
    'so_sat3': EarthlySite(-_debabyl(67,47,16.53), -_debabyl(22,57,35.97), 5188.,
                           typical_weather=weather_factory('toco')),
}
# Take LAT as default.  It makes most sense to use a single position
# for all telescopes and absorb the discrepancy into the pointing
# model as a few seconds of base tilt.
SITES['so'] = SITES['so_lat']
SITES['_default'] = SITES['so']
DEFAULT_SITE = 'so'  # deprecated, use SITES['_default']

# These definitions are used for the naive horizon -> celestial
# conversion.
ERA_EPOCH = 946684800 + 3600 * 12  # Noon on Jan 1 2000.
ERA_POLY = np.array([6.300387486754831,
                     4.894961212823756])  # Operates on "days since ERA_EPOCH".


class CelestialSightLine:
    """Carries a vector of celestial pointing data.

    Instantiation should result in the pointing quaternion being
    stored in self.Q.

    """
    @staticmethod
    def decode_site(site=None):
        """Convert site argument to a Site object.  The argument must be one
        of:

        - an EarthlySite object
        - a string, corresponding to the name of an internally known
          site
        - None, in which casethe default site info will be loaded.

        """
        if site is None:
            site = DEFAULT_SITE
        if isinstance(site, EarthlySite):
            return site
        if site in SITES:
            return SITES[site]
        raise ValueError("Could not decode %s as a Site." % site)

    @classmethod
    def naive_az_el(cls, t, az, el, roll=0., site=None, weather=None):
        """Construct a SightLine from horizon coordinates, az and el (in
        radians) and time t (unix timestamp).

        This will be off by several arcminutes... but less than a
        degree.  The weather is ignored.

        """
        site = cls.decode_site(site)
        assert isinstance(site, EarthlySite)

        self = cls()

        J = (t - ERA_EPOCH) / 86400
        era = np.polyval(ERA_POLY, J)
        lst = era + site.lon * DEG

        self.Q = (
            quat.euler(2, lst) *
            quat.euler(1, np.pi/2 - site.lat * DEG) *
            quat.euler(2, np.pi) *
            quat.euler(2, -az) *
            quat.euler(1, np.pi/2 - el) *
            quat.euler(2, roll)
            )
        return self

    @classmethod
    def az_el(cls, t, az, el, roll=None, site=None, weather=None, **kwargs):
        """Construct a SightLine from horizon coordinates.  This uses
        high-precision pointing. kwargs are passed to qpoint.Qpoint.

        """
        import qpoint  # https://github.com/arahlin/qpoint

        site = cls.decode_site(site)
        assert isinstance(site, EarthlySite)

        if isinstance(weather, str):
            if weather == 'typical':
                weather = site.typical_weather
            else:
                weather = weather_factory(weather)

        if weather is None:
            raise ValueError('High-precision pointing requires a weather '
                             'object.  Try passing \'toco\', \'typical\', or '
                             '\'vacuum\'.')

        self = cls()
        # Note that passing num_threads explicitly to qpoint will
        # cause openmp_set_thread_count to be called!
        qp = qpoint.QPoint(accuracy='high', fast_math=True, mean_aber=True,
                           rate_ref='always', **weather.to_qpoint(), **kwargs)

        az, el, t = map(np.asarray, [az, el, t])
        Q_arr = qp.azel2bore(az / DEG, el / DEG, None, None,
                             lon=site.lon, lat=site.lat, ctime=t)
        self.Q = quat.G3VectorQuat(Q_arr)

        # Apply boresight roll manually (note qpoint pitch/roll refer
        # to the platform rather than to something about the boresight
        # axis).  Regardless we need a pi rotation because of our
        # definition of boresight coordinates.
        if roll is None:
            self.Q *= quat.euler(2, np.pi)
        if roll is not None:
            self.Q *= quat.euler(2, np.pi + roll)

        return self

    @classmethod
    def for_lonlat(cls, lon, lat, psi=0):
        """Construct a SightLine directly from lonlat angles representing
        celestial coordinates.

        This is appropriate if you want a SightLine that takes focal
        plane coordinates directly to some celestial pointing.  In
        that case, lon and lat are RA and dec, and psi is the
        parallactic rotation of the focal plane on the sky at the
        target position.  psi=0 will correspond to focal plane "up"
        (+eta axis) mapping parallel to lines of longitude; psi > 0 is
        a clockwise rotation of the focal plane on the sky
        (i.e. opposite IAU Position Angle direction).

        """
        self = cls()
        self.Q = quat.euler(2, lon) * quat.euler(1, np.pi/2 - lat) * quat.euler(2, psi)
        return self

    @classmethod
    def for_horizon(cls, t, az, el, roll=None, site=None, weather=None):
        """Construct the trivial SightLine, where "celestial" coordinates are
        taken to simply be local horizon coordinates (without any
        effects of atmosphere).  Although this accepts site and
        weather arguments, they are ignored.

        """
        self = cls()
        az, el, t = map(np.asarray, [az, el, t])
        if roll is None:
            roll = 0.
        else:
            roll = np.asarray(roll)
        self.Q = (
            quat.euler(2, -az) *
            quat.euler(1, np.pi/2 - el) *
            quat.euler(2, roll)
            )
        return self

    def coords(self, fplane=None, output=None):
        """Get the celestial coordinates of each detector at each time.

        Arguments:
          fplane: A FocalPlane object representing the detector
            offsets and responses, or None
          output: An optional structure for holding the results.  For
            that to work, each element of output must support the
            buffer protocol.

        Returns:
          If fplane is not None, then the result will be
          [n_samp,{lon,lat,cos2psi,sin2psi}]. Otherwise it will
          be [n_det,n_Samp,{lon,lat,cos2psi,sin2psi}]
        """
        # Get a projector, in CAR.
        p = so3g.ProjEng_CAR_TQU_NonTiled((1, 1, 1., 1., 1., 1.))
        # Pre-process the offsets
        collapse = (fplane is None)
        if collapse:
            fplane = FocalPlane()
            if output is not None:
                output = output[None]
        output = p.coords(self.Q, fplane.quats, output)
        if collapse:
            output = output[0]
        return output

class FocalPlane:
    """This class represents the detector positions and intensity and
    polarization responses in the focal plane.

    This used to be a subclass of OrderedDict, but this was hard to
    generalize to per-detector polarization efficiency.
    """
    def __init__(self, quats=None, resps=None):
        if quats is None: quats = np.array([[1.0,0.0,0.0,0.0]])
        # Building them this way ensures that
        # quats will be an quat coeff array-2 and resps will be a numpy
        # array with the right shape, so we don't need to check
        # for this when we use FocalPlane later
        self.quats = np.array(quat.G3VectorQuat(quats))
        self.resps = np.ones((len(self.quats),2),np.float32)
        if resps is not None:
            self.resps[:] = resps
        if np.any(~np.isfinite(self.quats)):
            raise ValueError("nan/inf values in detector quaternions")
        if np.any(~np.isfinite(self.resps)):
            raise ValueError("nan/inf values in detector responses")
    @classmethod
    def from_xieta(cls, xi, eta, gamma=0, T=1, P=1, Q=1, U=0, hwp=False):
        # The underlying code wants polangle gamma and the T and P
        # response, but we support speifying these as the T, Q and U
        # response too. Writing it like this handles both cases, and
        # as a bonus(?) any combination of them
        gamma = gamma + np.arctan2(U,Q)/2
        P     = P * (Q**2+U**2)**0.5
        if hwp: gamma = -gamma
        # Broadcast everything to 1d
        xi, eta, gamma, T, P, _ = np.broadcast_arrays(xi, eta, gamma, T, P, [0])
        quats = quat.rotation_xieta(xi, eta, gamma)
        resps = np.ones((len(quats),2))
        resps[:,0] = T
        resps[:,1] = P
        return cls(quats, resps=resps)
    def __repr__(self):
        return "FocalPlane(quats=%s,resps=%s)" % (repr(self.quats), repr(self.resps))
    @property
    def ndet(self): return len(self.quats)
    def __getitem__(self, sel):
        return FocalPlane(quats=self.quats[sel], resps=self.resps[sel])

class Assembly:
    """This class groups together boresight and detector offset
    information for use by the Projectionist.

    The basic implementation simply attachs some number of detectors
    rigidly to a boresight.  But this abstraction layer could
    facilitate more complex arrangements, eventually.

    """
    def __init__(self, collapse=False):
        self.collapse = collapse

    @classmethod
    def attach(cls, sight_line, fplane):
        """Create an Assembly based on a CelestialSightLine and a FocalPlane.

        Args:
            sight_line (CelestialSightLine): The position and
                orientation of the boresight.  This can just be a
                G3VectorQuat if you don't have a whole
                CelestialSightLine handy.
            fplane (FocalPlane): The "offsets" of each detector
                relative to the boresight, and their response to
                intensity and polarization
        """
        self = cls()
        if isinstance(sight_line, quat.G3VectorQuat):
            self.Q = sight_line
        else:
            self.Q = sight_line.Q
        self.fplane = fplane
        return self

    @classmethod
    def for_boresight(cls, sight_line):
        """Returns an Assembly where a single detector serves as a dummy for
        the boresight.

        """
        self = cls(collapse=True)
        self.Q = sight_line.Q
        self.fplane = FocalPlane()
        return self

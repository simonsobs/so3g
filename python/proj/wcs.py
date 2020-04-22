import so3g
from . import quat

import numpy as np

from .ranges import Ranges, RangesMatrix

class Projectionist:
    """This class assists with analyzing WCS information to populate data
    structures needed for accelerated pointing routines.

    The optimized pointing code only needs to know about the Native
    Celestial Coordinate system, because the projection plane and thus
    the pixel indices are defined relative to that.

    As in pixell, the code below uses the term "geometry" to refer to
    the combination of an astropy.WCS and object and a 2-d array
    shape.

    """
    @staticmethod
    def get_q(wcs):
        """Analyze a wcs object to compute the quaternion rotation from
        celestial to native spherical coordinates.

        """
        alpha0, delta0 = wcs.wcs.crval  # In degrees.
        if (wcs.wcs.phi0 == 0. and wcs.wcs.theta0 == 0.):
            # This is typical for cylindrical projections.
            assert((delta0 >= 0 and wcs.wcs.lonpole == 180.0) or
                   (delta0 <= 0 and wcs.wcs.lonpole ==   0.0))
            Q = (quat.euler(1,  delta0 * quat.DEG) *
                 quat.euler(2, -alpha0 * quat.DEG))
        elif (wcs.wcs.phi0 == 0. and wcs.wcs.theta0 == 90.):
            # This is typical for zenithal projections.
            assert(wcs.wcs.lonpole == 180.0)
            Q = (quat.euler(2, np.pi) *
                 quat.euler(1, (delta0 - 90)*quat.DEG) *
                 quat.euler(2, -alpha0 * quat.DEG))
        else:
            raise ValueError(f'Unimplemented NSC reference (phi0,theta0)='
                             '({wcs.wcs.phi0:.2f},{wcs.wcs.theta0:.2f})')
        return Q

    def __init__(self):
        self._q0 = None
        self.naxis = np.array([0, 0])
        self.cdelt = np.array([0., 0.])
        self.crpix = np.array([0., 0.])

    @classmethod
    def for_geom(cls, shape, wcs):
        self = cls()
        ax1, ax2 = wcs.wcs.lng, wcs.wcs.lat
        self.ndim = len(shape)
        # The axes are numbered from outside in...
        self.naxis = np.array([shape[self.ndim - ax1 - 1],
                               shape[self.ndim - ax2 - 1]], dtype=int)
        # Get just the celestial part.
        wcs = wcs.celestial

        # Extract the projection name (e.g. CAR)
        proj = [c[-3:] for c in wcs.wcs.ctype]
        assert(proj[0] == proj[1])
        proj_name = proj[0]  # Projection name
        self.proj_name = proj_name

        # Store the rotation to native spherical coordinates.
        self.q_celestial_to_native = self.get_q(wcs)

        # Store the grid info.
        self.cdelt = np.array(wcs.wcs.cdelt) * quat.DEG
        self.crpix = np.array(wcs.wcs.crpix)

        return self

    @classmethod
    def for_map(cls, emap, wcs=None):
        if wcs is None:
            wcs = emap.wcs
        return cls.for_geom(emap.shape, wcs)

    @classmethod
    def for_source_at(cls, alpha0, delta0, gamma0=0.,
                      proj_name='TAN'):
        """Return a projection that places the specified source position at
        the North Pole."""

        self = cls()
        self.proj_name = proj_name
        assert(gamma0 == 0.)
        self.q_celestial_to_native = (
            quat.euler(2, np.pi)
            * quat.euler(1, (delta0 - 90)*quat.DEG)
            * quat.euler(2, -alpha0 * quat.DEG))
        return self

    def get_pixelizor(self):

        """Returns the so3g.Pixelizor appropriate for use with the configured
        geometry.

        """
        # All these casts are required because boost-python doesn't
        # like numpy scalars.
        return so3g.Pixelizor2_Flat(int(self.naxis[1]), int(self.naxis[0]),
                                    float(self.cdelt[1]), float(self.cdelt[0]),
                                    float(self.crpix[1]), float(self.crpix[0]))

    def get_ProjEng(self, comps='TQU', proj_name=None, get=True,
                    instance=True):
        """Returns an so3g.ProjEng object appropriate for use with the
        configured geometry.

        """
        if proj_name is None:
            proj_name = self.proj_name
        projeng_name = f'ProjEng_{proj_name}_{comps}'
        if not get:
            return projeng_name
        try:
            projeng_cls = getattr(so3g, projeng_name)
        except AttributeError:
            raise ValueError(f'There is no projector implemented for '
                             'pixelization "{proj_name}", components '
                             '"{comps}" (tried "{projeng_name}").')
        if not instance:
            return projeng_cls
        return projeng_cls(self.get_pixelizor())

    def _get_cached_q(self, new_q0):
        if new_q0 is not self._q0:
            self._q0 = new_q0
            self._qv = self.q_celestial_to_native * self._q0
        return self._qv

    def _guess_comps(self, map_shape):
        if len(map_shape) != 3:
            raise ValueError('Cannot guess components based on '
                             'map with %i!=3 dimensions!' % len(map_shape))
        if map_shape[0] == 1:
            return 'T'
        elif map_shape[0] == 2:
            return 'QU'
        elif map_shape[0] == 3:
            return 'TQU'
        raise ValueError('No default components for ncomp=%i; '
                         'set comps=... explicitly.' % map_shape[0])

    def get_pixels(self, assembly):
        """Get the pixel indices for the provided pointing Assembly.  For each
        detector, an int32 array of shape [n_time] is returned.

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('TQU')
        q1 = self._get_cached_q(assembly.Q)
        return projeng.pixels(q1, assembly.dets, None)

    def get_pointing_matrix(self, assembly):
        """Get the pointing matrix information, which is to say both the pixel
        indices and the "spin projection factors" for the provided
        pointing Assembly.  Returns (pixel_indices, spin_proj) where
        pixel_indices is a list of int32 arrays of shape [n_time] and
        spin_proj is a list of float arrays of shape [n_time,
        n_component].  As an alternative to on-the-fly computation,
        this information can be cached and used for very fast
        projection operations.

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('TQU')
        q1 = self._get_cached_q(assembly.Q)
        return projeng.pointing_matrix(q1, assembly.dets, None, None)

    def get_coords(self, assembly, use_native=False):
        """Get the spherical coordinates for the provided pointing Assembly.
        For each detector, a float64 array of shape [n_time,4] is
        returned.  In the right-most index, the first two components
        are the longitude and latitude in radians.  The next two
        components are the cosine and sine of the parallactic angle.

        This routine uses a trivial CAR projection to return results
        from the celestial coordinate system (i.e. (lon,lat) =
        (RA,dec)), and the parallactic angle will be relative to that
        projection.  If you are interested in the parallactic angle of
        the native spherical coordinate system (e.g. if you're doing a
        tangent plane projection), make a second call specifying
        use_native=True.  In this case you might also take a look at
        the get_planar() routine.

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('TQU', 'CAR')
        if use_native:
            q1 = self._get_cached_q(assembly.Q)
        else:
            q1 = assembly.Q
        return projeng.coords(q1, assembly.dets, None)

    def get_planar(self, assembly):
        """Get projection plane coordinates for all detectors at all times.
        For each detector, a float64 array of shape [n_time,4] is
        returned.  The first two elements are the x and y projection
        plane coordiantes, similar to the "intermediate world
        coordinates", in FITS language.  Insofar as FITS ICW has units
        of degrees, these coordinates have units of radians.  Indices
        2 and 3 carry the cosine and sine of the detector parallactic
        angle.

        """
        q1 = self._get_cached_q(assembly.Q)
        projeng = self.get_ProjEng('TQU')
        return projeng.coords(q1, assembly.dets, None)

    def get_prec_omp(self, assembly):
        """Perform a quick analysis of the pointing in order to enable OMP in
        tod-to-map operations.  Returns a special object that can be
        passed as the omp= argument of to_map and to_weight_map.
        (Other operations can use OMP without this precomputed object,
        because there are no thread-safety issues.)

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('T')
        q1 = self._get_cached_q(assembly.Q)
        omp_ivals = projeng.pixel_ranges(q1, assembly.dets)
        return RangesMatrix([RangesMatrix(x) for x in omp_ivals])

    def to_map(self, signal, assembly, dest_map=None, omp=None, comps=None):
        """Project signal into a map.

        Arguments:
          signal (Signal-like): The signal to project.
          dest_map (Map-like): The map into which to accumulate the
            projected signal.  If None, a map will be initialized internally.
          comps: The projection component string, e.g. 'T', 'QU',
            'TQU'.
          omp (ProjectionOmpData): The OMP information (returned by
            get_prec_omp), if OMP acceleration is to be used.

        See class documentation for description of standard arguments.

        """
        if dest_map is None and comps is None:
            raise ValueError("Provide an output map or specify component of "
                             "interest (e.g. comps='TQU').")
        if comps is None:
            comps = self._guess_comps(dest_map.shape)
        projeng = self.get_ProjEng(comps)
        q1 = self._get_cached_q(assembly.Q)
        if omp is None:
            map_out = projeng.to_map(
                dest_map, q1, assembly.dets, signal, None)
        else:
            map_out = projeng.to_map_omp(
                dest_map, q1, assembly.dets, signal, None, omp)
        return map_out

    def to_weights(self, assembly, dest_map=None, omp=None, comps=None):
        """Project pointing into a weights map.

        Arguments:
          dest_map (Map-like): The weights map into which to
            accumulate the projected signal.  If None, a map will be
            initialized internally.
          comps: The projection component string, e.g. 'T', 'QU',
            'TQU'.
          omp (ProjectionOmpData): The OMP information (returned by
            get_prec_omp), if OMP acceleration is to be used.

        See class documentation for description of standard arguments.

        """
        if dest_map is None and comps is None:
            raise ValueError("Provide an output map or specify component of "
                             "interest (e.g. comps='TQU').")
        if comps is None:
            assert(dest_map.shape[0] == dest_map.shape[1])
            comps = self._guess_comps(dest_map.shape[1:])
        projeng = self.get_ProjEng(comps)
        q1 = self._get_cached_q(assembly.Q)
        if omp is None:
            map_out = projeng.to_weight_map(
                dest_map, q1, assembly.dets, None, None)
        else:
            map_out = projeng.to_weight_map_omp(
                dest_map, q1, assembly.dets, None, None, omp)
        return map_out

    def from_map(self, src_map, assembly, signal=None, comps=None):
        """De-project from a map, returning a Signal-like object.

        Arguments:
          src_map (Map-like): The map from which to sample.
          signal (Signal-like): The object into which to accumulate
            the signal.  If not provided, a suitable object will be
            created and initialized to zero.
          comps: The projection component string, e.g. 'T', 'QU',
            'TQU'.
          omp: The OMP information (returned by get_prec_omp), if OMP
            acceleration is to be used.

        See class documentation for description of standard arguments.

        """
        if src_map.ndim == 2:
            # Promote to (1,...)
            src_map = src_map[None]
        if comps is None:
            comps = self._guess_comps(src_map.shape)
        projeng = self.get_ProjEng(comps)
        q1 = self._get_cached_q(assembly.Q)
        signal_out = projeng.from_map(
            src_map, q1, assembly.dets, signal, None)
        return signal_out

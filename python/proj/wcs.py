import so3g
from . import quat

import numpy as np

from .ranges import Ranges, RangesMatrix

class Projectionist:
    """This class assists with analyzing WCS information to populate data
    structures needed for accelerated pointing routines.

    On instantiation, it carries information about the relation
    between celestial spherical coordinates and the Native Spherical
    coordinates of the projection, as and also the pixelization scheme
    for the projection.

    As in pixell, the code and discussion here uses the term
    "geometry" to refer to the combination of an astropy.WCS object
    and a 2-d array shape.

    When coordinate computation or projection routines are called, a
    ProjEng (a wrapped C++ object) is instantiated to perform the
    accelerated computations.  The spin-composition (i.e. is this a
    'T' map or a 'TQU' map) must be specified at this time, to
    instantiate the correct accelerated object.

    Some method parameters are common to many functions and are documented
    here for consistency:

    * ``assembly`` - an Assembly, providing a focal plane (quaternion
      offsets for each detector) as well as a boresight vector
      (quaternion for each time sample).
    * ``comps`` - a string specifying the Stokes components of the
      map, for example 'T' or 'TQU'.  When absent, this will be
      guessed from the map shape; with 1|2|3 mapping to 'T'|'QU'|'TQU'
      respectively.
    * ``proj_name`` - a string specifying a projection.  The
      nomenclature is mostly the same as the FITS CTYPE identifiers.
      Accepted values: ARC, CAR, CEA, TAN, ZEA, Flat, Quat.
    * ``threads`` - a RangesMatrix with shape
      (n_threads,n_dets,n_samps), used to specify which samples should
      be treated by each thread in TOD-to-map operations.  Such
      objects should satisfy the condition that
      threads[x,j]*threads[y,j] is the empty Range; i.e. each
      detector-sample is assigned to exactly one thread.

    Attributes:
        naxis: 2-element integer array specifying the map shape (for
            the 2 celestial axes).
        cdelt: 2-element float array specifying the pixel pitch.
        crpix: 2-element float array specifying the pixel coordinates
            of the reference point.
        proj_name: string, name of the projection.
        q_celestial_to_native: quaternion rotation taking celestial
            coordinates to the native spherical coordinates of the
            projection.

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
        """Return a Projectionist for use with the specified "geometry".

        The shape and wcs are the core information required to prepare
        a Projectionist, so this method is likely to be called by
        other constructors.

        """
        self = cls()
        ax1, ax2 = wcs.wcs.lng, wcs.wcs.lat
        ndim = len(shape)
        # The axes are numbered from outside in...
        self.naxis = np.array([shape[ndim - ax1 - 1],
                               shape[ndim - ax2 - 1]], dtype=int)
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
        """Return a Projectionist for use with maps having the same geometry
        as the provided enmap.

        Args:
          emap: enmap from which to extract shape and wcs information.
            It is acceptable to pass a bare ndarray here (or anything
            with shape attribute), provided that wcs is provided
            separately.
          wcs: optional WCS object to use instead of emap.wcs.

        """
        if wcs is None:
            wcs = emap.wcs
        return cls.for_geom(emap.shape, wcs)

    @classmethod
    def for_source_at(cls, alpha0, delta0, gamma0=0.,
                      proj_name='TAN'):
        """Return a pointing-only Projectionist where some particular
        equatorial position will be put at the North Pole of the
        Native spherical coordinates.

        """
        self = cls()
        self.proj_name = proj_name
        assert(gamma0 == 0.)
        self.q_celestial_to_native = (
            quat.euler(2, np.pi)
            * quat.euler(1, (delta0 - 90)*quat.DEG)
            * quat.euler(2, -alpha0 * quat.DEG))
        return self

    def _get_pixelizor_args(self):
        """Returns a tuple of arguments that may be passed to the ProjEng
        constructor to define the pixelization.

        """
        # All these casts are required because boost-python doesn't
        # like numpy scalars.
        return (int(self.naxis[1]), int(self.naxis[0]),
                float(self.cdelt[1]), float(self.cdelt[0]),
                float(self.crpix[1]), float(self.crpix[0]))

    def get_ProjEng(self, comps='TQU', proj_name=None, get=True,
                    instance=True):
        """Returns an so3g.ProjEng object appropriate for use with the
        configured geometry.

        """
        if proj_name is None:
            proj_name = self.proj_name
        projeng_name = f'ProjEng_{proj_name}_{comps}_NonTiled'
        if not get:
            return projeng_name
        try:
            projeng_cls = getattr(so3g, projeng_name)
        except AttributeError:
            raise ValueError(f'There is no projector implemented for '
                             f'pixelization "{proj_name}", components '
                             f'"{comps}" (tried "{projeng_name}").')
        if not instance:
            return projeng_cls
        return projeng_cls(self._get_pixelizor_args())

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

    def to_map(self, signal, assembly, output=None, det_weights=None,
               threads=None, comps=None):
        """Project signal into a map.

        Arguments:
          signal (Signal-like): The signal to project.
          output (Map-like): The map into which to accumulate the
            projected signal.  If None, a map will be initialized internally.
          det_weights (shape n_det): Weights to use for each detector.

        See class documentation for description of standard arguments.

        """
        if output is None and comps is None:
            raise ValueError("Provide an output map or specify component of "
                             "interest (e.g. comps='TQU').")
        if comps is None:
            comps = self._guess_comps(output.shape)
        projeng = self.get_ProjEng(comps)
        q1 = self._get_cached_q(assembly.Q)
        map_out = projeng.to_map(
            output, q1, assembly.dets, signal, det_weights, threads)
        return map_out

    def to_weights(self, assembly, output=None, det_weights=None,
                   threads=None, comps=None):
        """Project pointing into a weights map.

        Arguments:
          output (Map-like): The weights map into which to accumulate
            the projected signal.  If None, a map will be initialized
            internally.
          det_weights (shape n_det): Weights to use for each detector.

        See class documentation for description of standard arguments.

        """
        if output is None and comps is None:
            raise ValueError("Provide an output map or specify component of "
                             "interest (e.g. comps='TQU').")
        if comps is None:
            assert(output.shape[0] == output.shape[1])
            comps = self._guess_comps(output.shape[1:])
        projeng = self.get_ProjEng(comps)
        q1 = self._get_cached_q(assembly.Q)
        map_out = projeng.to_weight_map(
            output, q1, assembly.dets, det_weights, threads)
        return map_out

    def from_map(self, src_map, assembly, signal=None, comps=None):
        """De-project from a map, returning a Signal-like object.

        Arguments:
          src_map (Map-like): The map from which to sample.
          signal (Signal-like): The object into which to accumulate
            the signal.  If not provided, a suitable object will be
            created and initialized to zero.

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
            src_map, q1, assembly.dets, signal)
        return signal_out

    def assign_threads(self, assembly):
        """Get a thread assignment RangesMatrix, using simple information
        returned by the underlying pixelization.

        The returned object can be passed to the ``threads`` argument
        of the projection methods in this object.

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('T')
        q1 = self._get_cached_q(assembly.Q)
        omp_ivals = projeng.pixel_ranges(q1, assembly.dets, None)
        return RangesMatrix([RangesMatrix(x) for x in omp_ivals])

    def assign_threads_from_map(self, assembly, tmap):
        """Assign threads based on a map.

        The returned object can be passed to the ``threads`` argument
        of the projection methods in this object.

        See class documentation for description of standard arguments.

        Args:
          tmap: Map structure with shape (1,m,n) where each pixel,
            when converted to integer, species the index of the thread
            to which that pixel should be assigned.

        """
        projeng = self.get_ProjEng('T')
        q1 = self._get_cached_q(assembly.Q)
        omp_ivals = projeng.pixel_ranges(q1, assembly.dets, tmap)
        return RangesMatrix([RangesMatrix(x) for x in omp_ivals])

import so3g
from . import quat

import numpy as np

from .ranges import Ranges, RangesMatrix
from . import mapthreads

# For coordinate systems we use the following abbreviations:
#
# - DC: Detector coordinates
# - FP: Focal plane coordinates
# - CS: Celestial spherical coordinates (likely equatorial)
# - NS: Native spherical coordinates of the map projection
#
# Saying that a quaternion rotation q_B<A "takes A to B", where A and
# B are coordinate systems, means that if a vector has coordinates v_A
# in system A then its coordinates in system B are:
#
#      v_B = q_B<A v_A q_B<A*
#
# The rotation taking detector coordinates to native spherical
# coordinates is the product
#
#     q_NS<DC = q_NS<CS q_CS<FP q_FP<DC
#
# In this module, quaternion rotations are written as:
# - Q_CS<FP: Assembly.Q (vector in samples)
# - Q_FP<DC: Assembly.dets (vector in dets)
# - Q_NS<CS: Projectionst.q_celestial_to_native
#
# The Projectionist caches a vector of Q_NS<FP, which is then simply
# called "q_native".  This is usually the thing to pass to C++ level.


#: Valid settings for "interpol".  First entry is the default.
INTERPOLS = ['nearest', 'bilinear']

class _ProjectionistBase:
    """This is a base class to assist with populating data
    structures needed for accelerated pointing routines.
    This class should not be used directly; call instead
    :class:`Projectionist` for rectangular pixelizations
    or :class:`ProjectionistHealpix` for HEALPix.

    On instantiation, it carries information about the relation
    between celestial spherical coordinates and the Native Spherical
    coordinates of the projection, and also the pixelization scheme
    for the projection.


    The following symbols are used in docstrings to represent shapes:

    * ``n_dets`` - Number of detectors in the focal plane.
    * ``n_samps`` - Number of samples in the signal space of each detector.
    * ``n_threads`` - Number of threads to use in parallelized computation.
    * ``n_mapidx`` - Number of indices required to specify a pixel in
      a map.  For simple RectPix maps this is probably 2, with the first
      index specifying row and the second index specifying column.
      But for tiled maps n_mapidx is 3, with the first axis specifying
      the tile number. For HEALPix, n_mapidx is 1 for untiled and 2
      for tiled maps.
    * ``n_comp`` - Number of map components, such as Stokes components
      or other spin components.  This takes the value 3, for example,
      if projection into TQU space has been requested.

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
    * ``threads`` - the thread assignment, consisting of a list of
      RangesMatrix objects.  Each RangesMatrix object must have shape
      (n_threads, n_dets, n_samps).  The n_threads does not need to
      be the same for every entry in the list.  In TOD-to-map
      operations, each entry of this list is processed fully before
      proceeding to the next one.  Each entry "ranges" is processed
      using (up to) the specified number of threads, such that thread
      i performs operations only on the samples included in
      ranges[i,:,:].  Most thread assignment routines in this module
      will return a list of two RangesMatrix objects,
      [ranges_parallel, ranges_serial].  The first item represents the
      part of the computation that can be done in parallel, and has
      shape (n_threads, n_dets, n_samps).  The ranges_serial object
      has shape (1, n_dets, n_samps) and represents any samples that
      need to be treated in a single thread.  The ranges_serial is
      only non-trivial when interpolation is active.

    """

    def __init__(self):
        raise NotImplementedError("Use child class Projectionist or ProjectionistHealpix")

    def _get_pixelizor_args(self):
        raise NotImplementedError("Use child class Projectionist or ProjectionistHealpix")

    def get_ProjEng(self, comps='TQU', proj_name=None, get=True,
                    instance=True, interpol=None):
        """Returns an so3g.ProjEng object appropriate for use with the
        configured geometry.

        """
        if proj_name is None: proj_name = self.proj_name
        tile_suffix = 'Tiled' if self.tiling else 'NonTiled'

        # Interpolation mode
        if interpol is None:
            interpol = self.interpol
        if interpol in ["nn", "nearest"]:
            interpol_suffix = ""
        elif interpol in ["lin", "bilinear"]:
            interpol_suffix = "_Bilinear"
        else:
            raise ValueError("ProjEng interpol '%s' not recognized" % str(interpol))

        projeng_name = f'ProjEng_{proj_name}_{comps}_{tile_suffix}{interpol_suffix}'
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

    def _cache_q_fp_to_native(self, q_fp_to_celestial):
        """Get q_fp_to_native for argument q_fp_to_celestial, and cache the
        result, and return that result later if called with the same
        argument.

        """
        if q_fp_to_celestial is not self._q_fp_to_celestial:
            self._q_fp_to_native = self.q_celestial_to_native * q_fp_to_celestial
            self._q_fp_to_celestial = q_fp_to_celestial
        return self._q_fp_to_native

    def _guess_comps(self, map_shape):
        if map_shape[0] == 1:
            return 'T'
        elif map_shape[0] == 2:
            return 'QU'
        elif map_shape[0] == 3:
            return 'TQU'
        raise ValueError('No default components for ncomp=%i; '
                         'set comps=... explicitly.' % map_shape[0])

    def _get_map_shape(self, imap):
        """
        Get map shape, supporting "bare" tiled maps (lists of None or tile array).
        For a list of tiles it returns the *tile shape*
        """
        if isinstance(imap, list): # Assume tile list
            for tile in imap:
                if tile is not None:
                    return tile.shape
            raise ValueError("map has no non-None entries")
        else:
            return imap.shape

    def get_pixels(self, assembly):
        """Get the pixel indices for the provided pointing Assembly.  For each
        detector, an int32 array of shape (n_samps, n_mapidx) is
        returned.  The value of the first slot of the second axis will
        be -1 if and only if the sample's projected position is
        outside the map footprint.

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('TQU')
        q_native = self._cache_q_fp_to_native(assembly.Q)
        return projeng.pixels(q_native, assembly.fplane.quats, None)

    def get_pointing_matrix(self, assembly):
        """Get the pointing matrix information, which is to say both the pixel
        indices and the "spin projection factors" for the provided
        pointing Assembly.  Returns (pixel_indices, spin_proj) where
        pixel_indices is as returned by get_pixels() and spin_proj is
        a list of float arrays of shape [n_samps, n_comp].  As an
        alternative to on-the-fly computation, this information can be
        cached and used for very fast projection operations.

        See class documentation for description of standard arguments.

        """
        projeng = self.get_ProjEng('TQU')
        q_native = self._cache_q_fp_to_native(assembly.Q)
        return projeng.pointing_matrix(q_native, assembly.fplane.quats, assembly.fplane.resps, None, None)

    def get_coords(self, assembly, use_native=False, output=None):
        """Get the spherical coordinates for the provided pointing Assembly.
        For each detector, a float64 array of shape [n_samps,4] is
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
            q_native = self._cache_q_fp_to_native(assembly.Q)
        else:
            q_native = assembly.Q
        return projeng.coords(q_native, assembly.fplane.quats, output)

    def get_planar(self, assembly, output=None):
        """Get projection plane coordinates for all detectors at all times.
        For each detector, a float64 array of shape [n_samps,4] is
        returned.  The first two elements are the x and y projection
        plane coordiantes, similar to the "intermediate world
        coordinates", in FITS language.  Insofar as FITS ICW has units
        of degrees, these coordinates have units of radians.  Indices
        2 and 3 carry the cosine and sine of the detector parallactic
        angle.

        """
        q_native = self._cache_q_fp_to_native(assembly.Q)
        projeng = self.get_ProjEng('TQU')
        return projeng.coords(q_native, assembly.fplane.quats, output)

    def zeros(self, super_shape):
        """Return a map, filled with zeros, with shape (super_shape,) +
        self.shape.  For tiled maps, this will be a list of map tiles
        (with shape (super_shape, ) + tile_shape.  If any tiles are
        not active, None is returned in the corresponding slots.

        """
        projeng = self.get_ProjEng('T')
        return projeng.zeros(super_shape)

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
            comps = self._guess_comps(self._get_map_shape(output))
        projeng = self.get_ProjEng(comps)
        q_native = self._cache_q_fp_to_native(assembly.Q)
        map_out = projeng.to_map(output, q_native, assembly.fplane.quats,
            assembly.fplane.resps, signal, det_weights, threads)
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
            shape = self._get_map_shape(output)
            assert(shape[0] == shape[1])
            comps = self._guess_comps(shape[1:])
        projeng = self.get_ProjEng(comps)
        q_native = self._cache_q_fp_to_native(assembly.Q)
        map_out = projeng.to_weight_map(output, q_native, assembly.fplane.quats,
            assembly.fplane.resps, det_weights, threads)
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
        if comps is None:
            comps = self._guess_comps(self._get_map_shape(src_map))
        projeng = self.get_ProjEng(comps)
        q_native = self._cache_q_fp_to_native(assembly.Q)
        signal_out = projeng.from_map(src_map, q_native, assembly.fplane.quats,
            assembly.fplane.resps, signal)
        return signal_out

    def assign_threads(self, assembly, method='domdir', n_threads=None):
        """Get a thread assignment RangesMatrix.  Different algorithms can be
        chosen, though this method does not provide fine-grained
        control of algorithm parameters and instead seeks to apply
        reasonable defaults.

        The methods exposed here are:

        - ``'simple'``: Divides the map into n_threads horizontal
          bands, and assigns the pixels in each band to a single
          thread.  This tends to be bad for scans that aren't
          horizontally correlated, or that don't have equal weight
          throughout the map.  The 'domdir' algorithm addresses both
          of those problems.
        - ``'domdir'``: For constant elevation scans, determines the
          dominant scan direction in the map and divides it into bands
          parallel to that scan direction.  The thickness of bands is
          adjusted so roughly equal numbers of samples fall into each
          band.
        - ``'tiles'``: In Tiled projections, assign each tile to a
          thread.  Some effort is made to balance the total number of
          samples over the threads.

        The returned object can be passed to the ``threads`` argument
        of the projection methods in this class.

        See class documentation for description of standard arguments.

        Args:
          method (str): Algorithm to use.

        """
        if method not in THREAD_ASSIGNMENT_METHODS:
            raise ValueError(f'No thread assignment method "{method}"; '
                             f'expected one of {THREAD_ASSIGNMENT_METHODS}')
        n_threads = mapthreads.get_num_threads(n_threads)

        if method == 'simple':
            projeng = self.get_ProjEng('T')
            q_native = self._cache_q_fp_to_native(assembly.Q)
            omp_ivals = projeng.pixel_ranges(q_native, assembly.fplane.quats, None, n_threads)
            return wrap_ivals(omp_ivals)

        elif method == 'domdir':
            fplane_rep = assembly.fplane[::100]
            if (self.tiling is not None) and (self.active_tiles is None):
                tile_info = self.get_active_tiles(assembly)
                active_tiles = tile_info['active_tiles']
                self.active_tiles = active_tiles
            else:
                active_tiles = None
            return mapthreads.get_threads_domdir(
                assembly.Q, assembly.fplane, shape=self.naxis[::-1], wcs=self.wcs,
                tile_shape=self.tile_shape, active_tiles=active_tiles,
                n_threads=n_threads, fplane_rep=fplane_rep)

        elif method == 'tiles':
            tile_info = self.get_active_tiles(assembly, assign=n_threads)
            self.active_tiles = tile_info['active_tiles']
            return tile_info['group_ranges']

        raise RuntimeError(f"Unimplemented method: {method}.")

    def assign_threads_from_map(self, assembly, tmap, n_threads=None):
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
        q_native = self._cache_q_fp_to_native(assembly.Q)
        n_threads = mapthreads.get_num_threads(n_threads)
        omp_ivals = projeng.pixel_ranges(q_native, assembly.fplane.quats, tmap, n_threads)
        return wrap_ivals(omp_ivals)

    def get_active_tiles(self, assembly, assign=False):
        """For a tiled Projection, figure out what tiles are hit by an
        assembly.

        See class documentation for description of standard arguments.

        Args:
          assign (bool or int): If True or an int, then the function
            will also compute a thread assignment, based on assigning
            each tile to a particular thread.  If this is an int, it
            sets the thread count; if it is simply True then
            OMP_NUM_THREADS is used.

        Returns:
          dict : The basic analysis results are:

            - ``'active_tiles'`` (list of int): sorted, non-repeating
              list of tiles that are hit by the assembly.
            - ``'hit_counts'`` (list of int): the number of hits in
              each tile returned in 'active_tiles', respectively.

          If the thread assignment took place then the dict will also
          contain:

            - ``'group_tiles'`` (list of lists of ints): There is one
              entry per thread, and the entry lists the tiles that
              have been assigned to that thread.
            - ``'group_ranges'`` (RangesMatrix): The thread
              assignments; shape (n_thread, n_det, n_samp).

        """
        if self.tiling is None:
            raise RuntimeError("This Projectionist not set up for Tiled maps.")
        projeng = self.get_ProjEng('T')
        q_native = self._cache_q_fp_to_native(assembly.Q)
        # This returns a G3VectorInt of length n_tiles giving count of hits per tile.
        hits = np.array(projeng.tile_hits(q_native, assembly.fplane.quats))
        tiles = np.nonzero(hits)[0]
        hits = hits[tiles]
        if assign is True:
            assign = so3g.useful_info()['omp_num_threads']
        if assign > 0:
            group_n = np.array([0 for g in range(assign)])
            group_tiles = [[] for _ in group_n]
            cands = sorted(list(zip(hits, tiles)), reverse=True) # [(1000, 12), (900, 4), ...]
            for _n, _t in cands:
                idx = group_n.argmin()
                group_n[idx] += _n
                group_tiles[idx].append(_t)
            imax = np.argmax(group_n)
            # max_ratio = group_n[imax] / np.mean(np.concatenate([group_n[:imax], group_n[imax+1:]]))
            # if len(group_n)>1 and max_ratio > 1.1:
            #     print(f"Warning: Threads poorly balanced. Max/mean hits = {max_ratio}")

            # Now paint them into Ranges.
            R = projeng.tile_ranges(q_native, assembly.fplane.quats, group_tiles)
            R = wrap_ivals(R)
            return {
                'active_tiles': list(tiles),
                'hit_counts': list(hits),
                'group_ranges': R,
                'group_tiles': group_tiles,
            }
        return {
            'active_tiles': list(tiles),
            'hit_counts': list(hits),
        }

    _ivals_format = 2

class Projectionist(_ProjectionistBase):
    """This class assists with analyzing WCS information to populate data
    structures needed for accelerated pointing routines for rectangular
    pixelization.

    See also the methods and parameters defined in the base class
    :class:`_ProjectionistBase`.

    As in pixell, the code and discussion here uses the term
    "geometry" to refer to the combination of an astropy.WCS object
    and a 2-d array shape.

    Some common method parameters specific to rectangular pixelization
    are documented here for consistency.

    * ``proj_name`` - a string specifying a projection.  The
      nomenclature is mostly the same as the FITS CTYPE identifiers.
      Accepted values: ARC, CAR, CEA, TAN, ZEA, Flat, Quat.
    * ``shape`` - the shape describing the celestial axes of the map.
      (For tiled representation this specifies the parent footprint.)
    * ``wcs`` - the WCS describing the celestial axes of the map.
      Together with ``shape`` this is a geometry; see pixell.enmap
      documentation.
    * ``interpol``: How positions that fall between pixel centers will
      be handled. Options are "nearest" (default): Use Nearest
      Neighbor interpolation, so a sample takes the value of
      whatever pixel is closest; or "bilinear": linearly
      interpolate between the four closest pixels. bilinear is
      slower (around 60%) but avoids problems caused by a
      discontinuous model.

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
        tile_shape: 2-element integer array specifying the tile shape
            (this is the shape of one full-size sub-tile, not the
            decimation factor on each axis).
        tiling: a _Tiling object if the map is tiled, None otherwise.

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
        self._q_fp_to_native = None
        self._q_fp_to_celestial = None
        self.tile_shape = None
        self.active_tiles = None
        self.wcs = None
        self.proj_name = None
        self.q_celestial_to_native = None
        self.interpol = INTERPOLS[0]

        self.naxis = np.array([0, 0])
        self.cdelt = np.array([0., 0.])
        self.crpix = np.array([0., 0.])

    @property
    def tiling(self):
        if self.tile_shape is None:
            return None
        return _Tiling(self.naxis[::-1], self.tile_shape)

    @classmethod
    def for_geom(cls, shape, wcs, interpol=None):
        """Construct a Projectionist for use with the specified "geometry".

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
        self.wcs = wcs

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

        # Pixel interpolation mode
        if interpol is None:
            interpol = INTERPOLS[0]
        self.interpol = interpol

        return self

    @classmethod
    def for_map(cls, emap, wcs=None, interpol=None):
        """Construct a Projectionist for use with maps having the same
        geometry as the provided enmap.

        Args:
          emap: enmap from which to extract shape and wcs information.
            It is acceptable to pass a bare ndarray here (or anything
            with shape attribute), provided that wcs is provided
            separately.
          wcs: optional WCS object to use instead of emap.wcs.
        """
        if wcs is None:
            wcs = emap.wcs
        return cls.for_geom(emap.shape, wcs, interpol=interpol)

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

    @classmethod
    def for_tiled(cls, shape, wcs, tile_shape, active_tiles=True, interpol=None):
        """Construct a Projectionist for use with the specified geometry
        (shape, wcs), cut into tiles of shape tile_shape.

        See class documentation for description of standard arguments.

        Args:
            tile_shape: tuple of ints, giving the shape of each tile.
            active_tiles: bool or list of ints.  Specifies which tiles
                should be considered active.  If True, all tiles are
                populated.  If None or False, this will remain
                uninitialized and attempts to generate blank maps
                (such as calls to zeros or to projection functions
                without a target map set) will fail.  Otherwise this
                must be a list of integers specifying which tiles to
                populate on such calls.

        """
        self = cls.for_geom(shape, wcs, interpol=interpol)
        self.tile_shape = np.array(tile_shape, 'int')
        if active_tiles is True:
            self.active_tiles = list(range(self.tiling.tile_count))
        elif active_tiles in [False, None]:
            self.active_tiles = None
        else:
            # Presumably a list of tile numbers.
            self.active_tiles = active_tiles
        return self

    def from_map(self, src_map, assembly, signal=None, comps=None):
        """De-project from a map, returning a Signal-like object.
        See parent class documentation for full description.

        """
        if src_map.ndim == 2:
            # Promote to (1,...)
            src_map = src_map[None]
        return super().from_map(src_map, assembly, signal, comps)

    def _get_pixelizor_args(self):
        """Returns a tuple of arguments that may be passed to the ProjEng
        constructor to define the pixelization.

        """
        # All these casts are required because boost-python doesn't
        # like numpy scalars.
        args = (int(self.naxis[1]), int(self.naxis[0]),
                float(self.cdelt[1]), float(self.cdelt[0]),
                float(self.crpix[1]), float(self.crpix[0]))
        if self.tiling:
            args += tuple(map(int, self.tile_shape))
            if self.active_tiles is not None:
                args += (self.active_tiles,)
        return args

    def _guess_comps(self, map_shape):
        if len(map_shape) != 3:
            raise ValueError('Cannot guess components based on '
                             'map with %i!=3 dimensions!' % len(map_shape))
        return super()._guess_comps(map_shape)


class ProjectionistHealpix(_ProjectionistBase):
    """Projectionist for Healpix maps.
    See base class :class:`_ProjectionistBase` for more methods
    and explanation of common method parameters.

    Attributes:
        nside: int, nside of the map, power of 2; 0 < nside <= 8192.
        nside_tile: None, int, or 'auto', nside of tiling.
                    Ntile will be 12*nside_tile**2. None for untiled.
                    If 'auto', an appropriate nside_tile will be computed
                    and set on calling :func:`compute_nside_tile`.
        ordering: str, 'NEST' or 'RING'. Only NEST supported for tiled maps.

    """
    def __init__(self):
        self._q_fp_to_native = None
        self._q_fp_to_celestial = None
        self.active_tiles = None
        self.proj_name = None
        self.q_celestial_to_native = quat.quat(1,0,0,0)
        self.interpol = 'nearest'
        self.tiling = None

        self.nside = None
        self.nside_tile = None
        self.ordering = 'NEST'  # 'RING' or 'NEST'. Only NEST can be used with tiles

    @classmethod
    def for_healpix(cls, nside, nside_tile=None, active_tiles=None, ordering='NEST', interpol=None):
        """Construct a Projectionist for Healpix maps.

        See class documentation for description of standard arguments.

        """
        self=cls()
        self.proj_name = 'HP'
        self.nside = nside
        self.active_tiles = active_tiles
        if interpol is not None and (interpol not in ['nearest', 'nn']):
            raise NotImplementedError("Only 'nearest' interpolation is supported for Healpix")
        self.interpol = 'nearest'

        self.ordering = ordering
        if ordering not in ['NEST', 'RING']:
            raise ValueError("ordering {ordering} should be 'NEST' or 'RING'")

        if nside_tile is not None:
            self.nside_tile = nside_tile
            self.tiling = True

        if ordering == 'RING' and self.nside_tile is not None:
            raise NotImplementedError("'RING' not supported for tiled maps")

        return self

    def compute_nside_tile(self, assembly, nActivePerThread=5, nThreads=None):
        """Automatically compute and set self.nside_tile for good balancing over threads.

        Arguments:
          nActivePerThread: int, how many active pixels do you want per thread (minimum)
          nThreads: int, number of threads to optimize for (takes from OMP_NUM_THREADS if None)

        """
        if self.nside_tile == 'auto':
            ## Estimate fsky
            nside_tile0 = min(4, self.nside)  # ntile = 192, for estimating fsky
            self.nside_tile = nside_tile0
            nActive = len(self.get_active_tiles(assembly)['active_tiles'])
            fsky = nActive / (12 * nside_tile0**2)
            if nThreads is None:
                nThreads = so3g.useful_info()['omp_num_threads']
            # nside_tile is smallest power of 2 satisfying nTile >= nActivePerThread * nthread / fsky
            self.nside_tile = int(2**np.ceil(0.5 * np.log2(nActivePerThread * nThreads / (12 * fsky))))
            self.nside_tile = min(self.nside_tile, self.nside)
        return self.nside_tile

    def get_active_tiles(self, assembly, assign=False):
        """For a tiled Projection, figure out what tiles are hit by an
        assembly. See parent class documentation for full description.

        """
        self.compute_nside_tile(assembly) # Set nside_tile if 'auto'
        return super().get_active_tiles(assembly, assign)

    def get_coords(self, assembly, use_native=False, output=None):
        """Get the spherical coordinates for the provided pointing Assembly.
        See parent class documentation for full description.

        """
        projeng = self.get_ProjEng('TQU')
        if use_native:
            q_native = self._cache_q_fp_to_native(assembly.Q)
        else:
            q_native = assembly.Q
        return projeng.coords(q_native, assembly.fplane.quats, output)

    def from_map(self, src_map, assembly, signal=None, comps=None):
        """De-project from a map, returning a Signal-like object.
        See parent class documentation for full description.

        """
        if self.nside_tile is None and src_map.ndim == 1:
            # Promote to (1,...)
            src_map = src_map[None]
        elif self.nside_tile is not None:
            if not isinstance(src_map, list):
                raise TypeError(f"Expected list for tiled map; src_map is {type(src_map)}")
            shape = self._get_map_shape(src_map)
            if len(shape) == 1:
                # Promote to (1,...)
                for itile in range(len(src_map)):
                    if src_map[itile] is not None:
                        src_map[itile] = src_map[itile][None]
        return super().from_map(src_map, assembly, signal, comps)

    def assign_threads(self, assembly, method=None, n_threads=None):
        """Get a thread assignment RangesMatrix.
        Available methods are ``'simple'`` and ``'tiles'``.
        See parent class documentation for full description.

        """
        if method is None:
            if self.nside_tile is None:
                method = 'simple'
            else:
                method = 'tiles'
        if (method not in THREAD_ASSIGNMENT_METHODS_HP) and \
           (method in THREAD_ASSIGNMENT_METHODS):
           raise ValueError(f'Thread assignment method "{method}" '
                            'not supported for ProjectionistHealpix. '
                            f'Expected one of {THREAD_ASSIGNMENT_METHODS_HP}.')
        return super().assign_threads(assembly, method, n_threads)

    def _get_pixelizor_args(self):
        """Returns a tuple of arguments that may be passed to the ProjEng
        constructor to define the pixelization.

        """
        if self.active_tiles is not None:
            active_tiles = list(map(int, self.active_tiles))
        else:
            active_tiles = None
        nside_tile = None
        if self.nside_tile is not None:
            nside_tile = int(self.nside_tile)

        args = (int(self.nside),
                int(self.ordering == 'NEST'),
                nside_tile,
                active_tiles)
        return args

    def _guess_comps(self, map_shape):
        if len(map_shape) != 2:
            raise ValueError('Cannot guess components based on '
                             'map with %i!=2 dimensions!' % len(map_shape))
        return super()._guess_comps(map_shape)


class _Tiling:
    """Utility class for computations involving tiled maps.

    """
    def __init__(self, shape, tile_shape):
        self.shape = shape[-2:]
        self.tile_shape = tile_shape
        nt0 = int(np.ceil(self.shape[0] / self.tile_shape[0]))
        nt1 = int(np.ceil(self.shape[1] / self.tile_shape[1]))
        self.super_shape = (nt0, nt1)
        self.tile_count = nt0 * nt1
    def tile_dims(self, tile):
        rowcol = self.tile_rowcol(tile)
        return (min(self.shape[0] - rowcol[0] * self.tile_shape[0], self.tile_shape[0]),
                min(self.shape[1] - rowcol[1] * self.tile_shape[1], self.tile_shape[1]))
    def tile_rowcol(self, tile):
        if tile >= self.tile_count:
            raise ValueError(f"Request for tile {tile} is out of bounds for "
                             f"super_shape={self.super_shape}")
        return (tile // self.super_shape[1], tile % self.super_shape[1])
    def tile_offset(self, tile):
        row, col = self.tile_rowcol(tile)
        return row * self.tile_shape[0], col * self.tile_shape[1]

def wrap_ivals(ivals):
    """Thread computation routines at C++ level return nested lists of
    Ranges objects; i.e. something like this::

      ivals = [
               [                             # thread assignments for first "bunch"
                 [Ranges, Ranges, ... ],     # for thread 0
                 [Ranges, Ranges, ... ],
                 ...
                 [Ranges, Ranges, ... ],     # for thread n-1.
               ],
               [                             # thread assignments for second "bunch"
                 [Ranges, Ranges, ... ],     # for thread 0
               ],
              ]

    This function wraps and returns each highest level entry into a
    RangesMatrix, i.e.::

       wrapped = [
               RangesMatrix(n_threads1, n_det, n_samp),
               RangesMatrix(n_threads2, n_det, n_samp),
       ]

    Currently all use cases have len(ivals) == 2 and n_threads2 = 1
    but the scheme is more general than that.

    """
    return [RangesMatrix([RangesMatrix(y) for y in x]) for x in ivals]

THREAD_ASSIGNMENT_METHODS = [
    'simple',
    'domdir',
    'tiles',
]

THREAD_ASSIGNMENT_METHODS_HP = [
    'simple',
    'tiles'
]

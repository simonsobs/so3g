====================================
Pointing and Projections (so3g.proj)
====================================

Overview
========

This module contains fast and/or flexible routines for pointing
computations that are required in the context of time-ordered data and
map-making.  The main focus is on rapid conversion from horizon
coordinates to celestial coordinates.  On-the-fly conversion to
projection plane coordinates and accumulation of detector data into
maps (without storing long position or pixel vectors) is also
supported.

The basic acceleration approach used here is to express the rotation
that takes a detector to the celestial sky as a product of a (fixed)
rotation that take each detector to some reference position in the
focal plane and a rotation that takes the reference position in the
focal plane to the celestial sky, as a function of time.  This is an
approximation insofar as non-linear effects such as atmospheric
refraction and aberration add an additional dependence on the horizon
and celestial pointing of a detector.  But having accounted for the
mean effects, at the reference position, the errors incurred on nearby
detectors in the focal plane may be small enough to ignore.

Dependencies
============

To take full advantage of this module, you will want to install these
packages:

- https://github.com/arahlin/qpoint : for high-precision astrometry.
- https://github.com/simonsobs/pixell : for rectilinear sky map
  operations.


Tutorial
========

Import
------

Because the Projection module has special requirements, it must be
explicitly imported for use.  It is likely you'll also need numpy::

  import so3g.proj
  import numpy as np

Site object
-----------

For anything involving both horizon and celestial coordinates, you
must define the telescope site.  You can construct an
:py:class:`so3g.proj.EarthlySite` with the longitude and latitude of
your site (see reference).  Or use the classmethod ``get_named`` to
retrieve some known position::

  # Get a Site object... how about the Atacama Cosmology Telescope
  site = so3g.proj.EarthlySite.get_named('act')

In addition to the site position on Earth, the object also contains
parameters for typical site weather, for refraction calculations.

Horizon to equatorial coordinates
---------------------------------

The :py:class:`so3g.proj.CelestialSightLine` class is used to store a
time-dependent rotation from focal plane coordinates to celestial
coordinates.  To create such an object based on the horizon
coordinates of a telescope boresight, you must pass in the az, el, and
time (as a unix timestamp), as well as a description of the telescope
site and the weather conditions::

  # Define DEG, carrying conversion from degrees to radians.
  DEG = np.pi / 180

  # Vectors of time, azimuth, and elevation.
  t = [1900000000.]
  az = [180. * DEG]
  el = [60. * DEG]

  # Construct a SightLine from this information.
  csl = so3g.proj.CelestialSightLine.az_el(t, az, el, site=site,
      weather='typical')

After instantiation, ``csl`` will have a member ``Q`` that holds the
rotation quaternion from focal plane to celestial coordinates.  The
``coords()`` method can decompose such quaternions into rotation
angles::

  >>> csl.Q
  spt3g.core.G3VectorQuat([(-0.0384775,0.941776,0.114177,0.313911)])

  >>> csl.coords()   
  array([[ 0.24261284, -0.9272257 , -0.99999913, -0.00131945]])
  
The ``coords()`` call returns an array with shape (n_time, 4); each
4-tuple contains values ``(lon, lat, cos(gamma), sin(gamma))``.  The
``lon`` and ``lat`` are the celestial longitude and latitude
coordinates (typically Right Ascension and Declination), in radians,
and gamma is the parallactic angle (the angle between the directions
of increasing declination and increasing elevation at that point, with
0 corresponding to North, increasing towards West).

You can get the vectors of RA and dec, in degrees like this::

  >>> ra, dec = csl.coords().transpose()[:2] / DEG
  >>> print(ra, dec)
  [13.90069189] [-53.12611929]


Pointing for many detectors
---------------------------

Create a FocalPlane object, with some detector positions and
orientations::

  names = ['a', 'b', 'c']
  x = np.array([-0.5, 0., 0.5]) * DEG
  y = np.zeros(3)
  gamma = np.array([0,30,60]) * DEG
  fp = so3g.proj.FocalPlane.from_xieta(names, x, y, gamma)


This particular function, ``from_xieta``, will apply the SO standard
coordinate definitions and store a rotation quaternion for each
detector.  FocalPlane is just a thinly wrapped OrderedDict, where the
detector name is the key and the value is the rotation quaternion::

  >>> fp['c']
  spt3g.core.quat(0.866017,0.00218168,0.00377878,0.499995)

At this point you could get the celestial coordinates for any one of
those detectors::

  # Get vector of quaternion pointings for detector 'a'
  q_total = csl.Q * fp['a']
  # Decompose q_total into lon, lat, roll angles
  ra, dec, gamma = so3g.proj.quat.decompose_lonlat(q_total)

As expected, these coordinates are close to the ones computed before,
for the boresight::

  >>> print(ra / DEG, dec / DEG)
  [13.90180432] [-53.6261252]

But the more expedient way to get pointing for multiple detectors is
to call ``coords()`` with the FocalPlane object as first argument::

  >>> csl.coords(fp)
  OrderedDict([('a', array([[ 0.24263226, -0.93595245, -0.99999911,
  -0.00133503]])), ('b', array([[ 0.24261284, -0.9272257 , -0.86536493,
  -0.50114224]])), ('c', array([[ 0.24259387, -0.91849895, -0.49887   ,
  -0.86667683]]))])

To be clear, ``coords()`` now returns a dictionary whose keys are the
detector names.  Each value is an array with shape (n_time,4), and at
each time step the 4 elements of the array are: ``(lon, lat,
cos(gamma), sin(gamma))``.


Projecting to Maps
------------------

Accelerated projection routines for various pixelizations have been
written using C++ templates.  Abstraction at the Python level is
provided by the :py:class:`so3g.proj.Projectionist` class, which
decodes the map WCS description in order to choose the right
accelerated pointing projection routines.

Let's start by creating a suitable map; we choose a supported
projection, and make sure it contains the particular points touched by
our example from above::

  from pixell import enmap
  shape, wcs = enmap.geometry([[-54*DEG, 16.*DEG], [-52*DEG, 12.*DEG]],
    res=.02*DEG, proj='car')

The ``shape`` is a simple tuple, and will be used like a numpy array
shape.  Because we typically expect array columns (fastest index) to
correspond to longitude, the ordering of the points is (dec, RA).  The
wcs is an astropy.wcs object, but decorated with a more useful repr
provided by pixell::

  >>> shape
  (150, 100)
  >>> wcs
  car:{cdelt:[-0.02,0.02],crval:[14,0],crpix:[101,2701]}

We can do lots of stuff with this ``(shape, wcs)`` pair... see the
`pixell documentation <https://pixell.readthedocs.io/>`__ .  One thing
we can do is get a map whose values provide the coordinates of each
pixel::

  pmap = enmap.posmap(shape, wcs)

You will see the shape is (2,100,200); the first dimension is because
there is a separate celestial map for the dec and the RA values.
Let's project those maps back into the time domain::

  p = so3g.proj.Projectionist.for_geom(shape, wcs)
  asm = so3g.proj.Assembly.attach(csl, fp)
  pix_dec = p.from_map(pmap[0], asm)
  pix_ra = p.from_map(pmap[1], asm)

Inspecting the values, we see they are close to the coordinates we
expect.  That makes sense, as these values correspond to the values
from the pixels in pmap, which are the coordinates of the centers of
the pixels::

  >>> [x/DEG for x in pix_ra]
  [array([13.88], dtype=float32), array([13.88], dtype=float32),
  array([13.88], dtype=float32)]
  >>> [x/DEG for x in pix_dec]
  [array([-53.6], dtype=float32), array([-53.100002], dtype=float32),
  array([-52.600002], dtype=float32)]

If you are not getting what you expect, you can grab the pixel indices
inferred by the projector -- perhaps your pointing is taking you off
the map (in which case the pixel indices would return value -1)::

  >>> p.get_pixels(asm)
  [array([4106], dtype=int32), array([9106], dtype=int32),
  array([14106], dtype=int32)]

Let's project signal into an intensity map::

  # Create dummy signal for our 3 detectors at 1 time points:
  signal = np.array([[1.], [10.], [100.]], dtype='float32')

  # Project into T-only map.
  map_out = p.to_map(signal, asm, comps='T')

Inspecting the map, we see our signal values occupy the three non-zero
pixels:

  >>> map_out.nonzero()
  (array([0, 0, 0]), array([20, 45, 70]), array([106, 106, 106]))
  >>> map_out[map_out!=0]
  array([  1.,  10., 100.])

If we run this projection again, but pass in this map as a starting
point, the signal will be added to the map:

  >>> p.to_map(signal, asm, dest_map=map_out, comps='T')
  array([[[0., 0., 0., ..., 0.]]])
  >>> map_out[map_out!=0]
  array([  2.,  20., 200.])

If we instead want to treat the signal as coming from
polarization-sensitive detectors, we can request components ``'TQU'``:

  map_pol = p.to_map(signal, asm, comps='TQU')

Now a 3-dimensional output map is created, with values binned
according to the projected detector angle on the sky::

  >>> map_pol.shape
  (3, 100, 200)
  >>> map_pol[:,45,106]
  array([10.        ,  4.97712898,  8.67341805])

For the most basic map-making, the other useful operation is the
``to_weights_map()`` method.  This is used to compute the weights
matrix in each pixel, returned as a (n_comp,n_comp,n_dec,n_ra) map.
Mathematically this corresponds to:

.. math::

   W = s^T\,P\,P^T\,s

where *s* is the time-ordered signal matrix with value 1 at every
position.  (Future extensions of this function will provide more
control over *s*.)

The weights map can be obtained like this::

  weight_out = p.to_weights(asm, comps='TQU')

The shape corresponds to a 3x3 matrix at each pixel; but notice only
the upper diagonal has been filled in, for efficiency reasons...::

  >>> weight_out.shape
  (3, 3, 100, 200)
  >>> weight_out[:,45,106]
  array([[1.        , 0.49771291, 0.86734182],
         [0.        , 0.24771814, 0.43168721],
         [0.        , 0.        , 0.75228184]])

OpenMP
------

The three routines of the Projectionist that move data between map and
time-domain are: ``from_map``, ``to_map``, and ``to_weights``.  The
``from_map`` function will automatically use OMP threading.  The
``to_map`` and ``to_weights`` functions need to be instructed on how
to use OpenMP safely, or else they will default to a single thread
implementation.

The reason that thread-safety is non-trivial for ``to_map`` and
``to_weights`` is that, depending on how the job is threaded, multiple
threads might try to update to the same pixel at the ~same time.  This
is a race condition that must be dealt with carefully.  Typical
approaches would be:

- Maintain separate copies of the output map for each thread, then add
  them together at the end.  (This is memory-intensive.)
- Use atomic update operations.  (This requires locking, which can be
  very slow to somewhat slow depending on the architecture.)

The approach that is currently implemented is to assign each pixel of
the map to a particular thread, and then for each detector to identify
all ranges of time samples that strike pixels belonging to each
thread.  Then each OpenMP thread loops over all detectors, but only
processes the ranges of samples that have been pre-identified as
belonging to that thread.  The result is that each thread touches a
disjoint set of pixels.

The precomputation needed to execute the above is non-trivial, and the
best way to assign pixels to threads depends on the particulars of the
scan pattern.  So OMP should be used with care.

The assignment of pixels to threads, and thus of sample-ranges to
threads, is encoded in a RangesMatrix object.  To get one, try
the ``Projectionist.get_prec_omp()`` method, then pass the result to
``to_map`` or ``from_map`` argument ``omp=``::

  omp_precomp = p.get_prec_omp(asm)
  map_pol2 = p.to_map(signal, asm, comps='TQU', omp=omp_precomp)

Inspecting::

  >>> omp_precomp
  RangesMatrix(4,3,1)
  >>> map_pol2[:,45,106]
  array([10.        ,  4.97712898,  8.67341805])


Class reference
===============

*The core classes from* ``so3g.proj`` *are auto-documented here.  If
you see a bunch of headings and no docstrings, then it's likely
because the Sphinx could not import so3g properly when building the
docs!*

Assembly
--------
.. autoclass:: so3g.proj.Assembly
   :members:

CelestialSightLine
------------------
.. autoclass:: so3g.proj.CelestialSightLine
   :members:

EarthlySite
-----------
.. autoclass:: so3g.proj.EarthlySite
   :members:

Projectionist
-------------
.. autoclass:: so3g.proj.Projectionist
   :members:

Ranges
------
.. autoclass:: so3g.proj.Ranges
   :members:

.. autoclass:: so3g.RangesInt32
   :members:

RangesMatrix
------------
.. autoclass:: so3g.proj.RangesMatrix
   :members:

Weather
-------
.. autoclass:: so3g.proj.Weather
   :members:

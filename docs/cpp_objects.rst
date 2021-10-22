C++ Objects
===========

Intervals
---------

The Intervals<T> objects have the following interface:

.. autoclass:: so3g.IntervalsDouble
   :members:

G3SuperTimestream
-----------------

G3SuperTimestream is for storing 2d arrays of data where the first
axis corresponds to named channels and the second axis indexes time.
The data arrays can have int32, int64, float32, or float64 data types.
It has configurable compression options in order to accomodate different
kinds of data.

For arrays of integers, lossless compression is enabled by default.
For arrays of floats, compression can be enabled that will be lossless
over a reduced dynamic and precision range.

Creating a G3SuperTimestream
````````````````````````````

When building a G3SuperTimestream, you must first populate the axis
information, then load in the data.  Here is an example::

  # imports
  from spt3g import core
  import so3g
  import numpy as np

  # The data we want to capture
  times = 1680000000 + 0.2 * np.arange(10000)
  names = ['a', 'b', 'c', 'd', 'e']
  data = (np.random.normal(size=(len(names), len(times))) * 256).astype('int32')

  # Creation of a G3SuperTimestream
  ts = so3g.G3SuperTimestream()
  ts.names = names
  ts.times = core.G3VectorTime(times * core.G3Units.s)
  ts.data = data


The object is now complete, and can be serialized.  In the default
configuration, arrays with int32 or int64 data types will be
compressed losslessly using a combination of FLAC and bzip.

Controlling Compression
```````````````````````

You can trigger compression of the data by calling ``.encode()``.  The
reference in ``.data`` is released, and a binary blob with compressed
data is saved internally.  If you try to access ``.data`` after
calling ``.encode()``, a new array will be created (by decompressing
the blob) and returned.

Compression will be triggered automatically on frame serialization
(i.e. when the frame is written to a file or network stream).  Because
serialization is a const operation, the binary blob of compressed data
is not stored in this case, and the reference to ``.data`` is not
released.  So there may be performance advantages to calling
``.encode()`` "manually" before passing your object through to
consumers that might want to use it in serialized form.

It is possible to tweak the compression algorithms, through the
``.options`` method, but this should be done with care.  For
compression evaluation and basic debugging one probably only wants to
use the highest level control, which simply enables or disables
compression::

  ts.options(enable=0)  # disable compression
  ts.options(enable=1)  # enable compression with default params

Two arguments allow some finer grain control over the FLAC and BZ2
algorithms and should not cause trouble (other than inefficiency) if
manipulated by the user:

  ``flac_level``
    The FLAC compression level, passed through to
    `FLAC__stream_encoder_set_compression_level`_.  Integer from 0 to
    8 with higher numbers corresponding to slower but potentially
    better compression.

  ``bz2_workFactor``
    The bzip2 workFactor, as described in `BZ2_bzCompressInit`_.  This
    has something to do with how soon the bz2 algorithm gives up on
    difficult (highly repetitive) data.

The additional arguments, `data_algo` and `times_algo`, are for
debugging and should not be messed with lightly.


.. _`FLAC__stream_encoder_set_compression_level`: https://xiph.org/flac/api/group__flac__stream__encoder.html#gae49cf32f5256cb47eecd33779493ac85
.. _`BZ2_bzCompressInit`: https://www.sourceware.org/bzip2/manual/manual.html#bzcompress-init

How to work with float arrays
`````````````````````````````

The G3SuperTimestream can also be used to carry non-integer data,
which is presented to the user as arrays of float32 or float64.  For
channels indexed by *i* and samples indexed by *t*, the non-integer
data array *x* will be discretized for compression and serialization
according to :math:`x_{it} = y_{it} \times q_i`, where *y* is an array
of integers and *q* are per-channel quanta (float64).

There are two ways to enter "float mode".  These are demonstrated in two examples.

**Example 1:** Populate a G3SuperTimestream with an integer array,
then apply a calibration factor (one per channel) using
``.calibrate``::

  # imports
  from spt3g import core
  import so3g
  import numpy as np

  # The data we want to capture
  times = 1680000000 + 0.2 * np.arange(10000)
  names = ['a', 'b', 'c', 'd', 'e']
  data = (np.random.normal(size=(len(names), len(times))) * 256).astype('int32')

  # Creation of a G3SuperTimestream
  ts = so3g.G3SuperTimestream()
  ts.names = names
  ts.times = core.G3VectorTime(times * core.G3Units.s)
  ts.data = data

  # Calibrate to, like, pW or something.
  pW_per_DAC = [1.23, 1.45, 1.89, 1.56, 1.01]
  ts.calibrate(pW_per_DAC)


**Example 2:** Populate a G3SuperTimestream by first assigning
discretization units (one per channel) to ``.quanta``, and then
assigning a float array to ``.data``::

  # imports
  from spt3g import core
  import so3g
  import numpy as np

  # The data we want to capture
  times = 1680000000 + 0.2 * np.arange(10000)
  names = ['a', 'b', 'c', 'd', 'e']
  data = (np.random.normal(size=(len(names), len(times))) * 256).astype('float64')

  # Creation of a G3SuperTimestream -- note we must set .quanta before
  # setting .data to a float array.
  ts = so3g.G3SuperTimestream()
  ts.names = names
  ts.times = core.G3VectorTime(times * core.G3Units.s)
  ts.quanta = 0.01 * np.ones(len(names))
  ts.data = data


This object is not suitable for lossless compression of arbitrary
float32 and float64 arrays.  The operator needs to have some idea of
what resolution must be preserved, i.e. how much rounding of arbitrary
float data is acceptable in the application.

Prior to compression the ``.data`` member will look like a float
array; but when compression is requested (or serialization occurs),
the following will happen:

- The array of floats *x* is converted to integers, *y = round(x /
  quantum)*.  If *x* is float32, then *y* will be packed into
  int32.  If *x* is float64, then *y* will be int64.
- The integer array *y* is compressed according to the lossless scheme
  used for integer data.

The precision of a particular non-zero float32 is approximately
:math:`2^{-23}` (1.2e-7) times its magnitude.  If we fix the precision
at 0.001, then the set of numbers that can be safely encoded by our
float32 scheme is all multiples of 0.001 between -8388.6 and +8388.6.

Equivalently for float64 the precision is :math:`2^{-52}` (2.2e-16)
times magnitude, so with a precision of 0.001 we have dynamic range of
about -4.5e12 to +4.5e12.

Working in C++
``````````````

If you need to construct many instances of G3SuperTimestream from
within C++, the method ``SetDataFromBuffer`` might be useful.  It will
copy data into a new numpy array from a C-ordered memory block,
allowing the caller to re-use the memory block.  A rough example is
presented below; see also the implementation of
``test_cxx_interface()`` in ``G3SuperTimestream.cxx``.

.. code-block:: c

   // Consider int32 array with 3 channels and 1000 samples.
   int shape[2] = {3, 1000};
   int typenum = NPY_INT32;

   // Use a flat buffer for storage.
   void *buf = calloc(shape[0] * shape[1], sizeof(int32_t));

   // ... fill up buf somehow ...

   // Create and manage a new G3SuperTimestream.
   auto ts = G3SuperTimestreamPtr(new G3SuperTimestream());

   // Set the channel names and timestamps.
   const char *chans[] = {"a", "b", "c"};
   ts->names = G3VectorString(chans, std::end(chans));
   ts->times = G3VectorTime();
   for (int i=0; i<n_samps; i++)
     ts->times.push_back(G3Time::Now());

   // Set compression options?
   // ts->Options(encode=0);

   // Set ts->data, by copying data from our buffer.
   ts->SetDataFromBuffer(buf, 2, shape, typenum, std::pair<int,int>(0, shape[1]));

   // Do something with ts...
   // writer->Process(ts);

   // Free what we allocated.
   free(buf);
  }


Interface autodoc
`````````````````

.. autoclass:: so3g.G3SuperTimestream
   :members:

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

It is possible to tweak the compression algorithm, but this should be
done with care.  By calling ``.options(data_algo=ALGO)`` you can set
the internal algorithm to one of the following:

- ALGO = 0: No compression, just store unmodified binary data.
- ALGO = 1: Use FLAC only.  This limits the dynamic range of the data
  to 24 bits and thus is not lossless if your input data exceeds this
  range.
- ALGO = 2: Use bzip only.  This is lossless but will not be efficient
  for "noisy" data.  There might be a use case here for arrays
  carrying slowly-changing bit-fields.
- ALGO = 3: Use FLAC+bzip (the default).  This is lossless and should
  perform well on noisy data.

On serialization, the ``.times`` vector is also compressed, using
bzip.  This can be disabled by passing ``.options(times_algo=0)``.


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
  import numpy as n

  # The data we want to capture
  times = 1680000000 + 0.2 * np.arange(10000)
  names = ['a', 'b', 'c', 'd', 'e']
  data = (np.random.normal(size=(len(names), len(times))) * 256).astype('float64')

  # Creation of a G3SuperTimestream -- note we must set .quanta before
  # setting .data to a float array.
  ts = so3g.G3SuperTimestream()
  ts.names = names
  ts.times = core.G3VectorTime(times * core.G3Units.s)
  ts.quanta = 0.01 * np.ones(len(chans))
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
   // ts->Options(data_algo=0);

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

Telescope Housekeeping Data (so3g.hk)
=====================================

The ``so3g.hk`` module contains code for writing and reading the
structured SO Housekeeping (HK) data format.  This format is carried
in spt3g G3 frames, but with a specific schema designed for SO use.
The detailed specification of the schema is contained in the
``tod2maps_docs/tod_format`` document.

If you're here you probably just want to read in some data, so we will
start with that.  Later on we go into details of the data model and
the interfaces for writing compliant HK files.

Reading HK Data
---------------

Typically one is working with an archive of HK files.  Loading the
data is a two-step process:

- Scan the files and cache meta-data about what fields are present at
  what times.
- Seek to the right parts of the right files and load the data for
  specific fields over a specific time range.

These two steps are accomplished with ``HKArchiveScanner`` and
``HKArchive``, respectively.  The rest of the section walks through an
example demonstrating basic usage of these two classes.

Scanning a set of files
```````````````````````

Here is how to instantiate an HKArchiveScanner and scan a bunch of
files::

  from so3g import hk
  import glob

  scanner = hk.HKArchiveScanner()

  files = glob.glob('/path/to/data/*.g3')
  for f in files:
      scanner.process_file(f)

  arc = scanner.finalize()

The returned object, ``arc``, is an instance of ``HKArchive``, and
that can be used to load the actual data vectors.

Reading the data
````````````````

To load data for a specific list of channels over a specific range of
times, call ``HKArchive.simple``::

    fields = ['LSA12.v_channel_1', 'LSA12.v_channel_2']
    data = arc.simple(fields)
    (t1, voltage1), (t2, voltage2) = data

At the end of this, ``t1`` and ``t2`` are vectors of timestamps and
``voltage1`` and ``voltage2`` contain the (ostensible) voltage
readings.

To restrict to some time range, pass unix timestamps as the next two
arguments (you can also pass them as keyword arguments, ``start=`` or
``end=``)::

    time_range = (1567800000, 1567900000)
    data = arc.simple(fields, time_range[0], time_range[1])

Abbreviating field names
````````````````````````

Full channel names in an observatory can be quite long
(e.g. "observatory.lat1_agg.LSA1234.Thermometer_29"), so this function
allows you to use shortened forms for the channel names in some
circumstances.  Suppose that the HK file set in the example above
contains only channel names beginning with ``'LSA12.'``.  Then the
function will understand if you ask for::

  fields = ['v_channel_1', 'v_channel_2']
  data = arc.simple(fields)

An error will be raised if the fields cannot be unambiguously matched
in the time range requested.  Note that the logic is only able to
include or exclude parts of the channel name between dots... so it
would not work to request ``fields = ['_1', '_2']``.  The caller can
suppress such matching by passing ``short_match=False``.

Co-sampling information
```````````````````````

Why is the method that gets data called ``simple``?  Because there is
a more sophisticated method called ``get_data`` that returns the data
in a somewhat more structured form.  This form is should be used when
you care about what fields are co-sampled.  See docstring.


Class references
````````````````

.. autoclass:: so3g.hk.HKArchive
   :members:

.. autoclass:: so3g.hk.HKArchiveScanner
   :members:

.. autoclass:: so3g.hk.HKScanner
   :members:

HK Data Types and File Structure [weak]
---------------------------------------

As of September 2019, all HK data that has been written is version 0
of the schema, which is only able to store vectors of double-precision
readings.  This will be extended substantially in version 1.

The HK file structures and versions are described in
https://github.com/simonsobs/tod2maps_docs/tod_format.

Writing HK Data [weak]
----------------------

Limited facilities exist to assist with creating valid HK data.  This
interface targets, especially, the OCS "aggregator" Agent.

.. autoclass:: so3g.hk.HKSessionHelper
   :members:


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

Loading Data Saved OCS/SO Style
---------------

One of the most basic things we might want to do is load data between
a time range. For `.g3` files that are saved by an OCS Aggregator, there
is a specific folder structure and file naming scheme. ``so3g.hk.load_range``
is writen to load data for a specified time frame saved in that style.

Example Use::

    from so3g.hk import load_range    
    data = load_range(start, stop, **kwargs)

Defining Time Ranges
```````````````````````

There are several options for defining start and stop. These options are passed
to ``so3g.hk.to_timestamp(some_time)`` to be parsed.

* ``datetime.datetime`` objects::

    import datetime as dt
    # tzinfo is likely needed if your computer is not in UTC
    start = dt.datetime(2020, 3, 12, 1, 23, tzinfo=dt.timezone.utc)
    
* Integers or Floats - these are assumed to be ctimes 

* Strings - assumed to be UTC dates and parsed by ``datetime.strptime``

    * '%Y-%m-%d'
    * '%Y-%m-%d %H:%M'
    * '%Y-%m-%d %H:%M:%S'
    * '%Y-%m-%d %H:%M:%S.%f'
    * Or submit your own with the ``str_format`` argument

Define Where to Look for Data
``````````````````````````````````````````````

* **Option 1** 

Use the ``data_dir = /path/to/ocs/hk/data`` keyword argument. This should be the
directory where all the first five digit ctime folders are located. ``load_range`` will 
look for the data there. 

* **Option 2** 

Set an environment variable ``export OCS_DATA_DIR = /path/to/ocs/hk/data``. 
``load_range`` will will automatically look there if it isn't over-ridden by option 1.

* **Option 3**

Use a configuration file. See Below.

Define Which Data to Load
``````````````````````````````````````````````

* **Option 1** 

No keyword arguments means ``load_range`` will return *every* field it can find. this
might take a very long time.

* **Option 2** 

Use the ``fields = [list, of, field, names]`` keyword argument. Example::
    
    fields = [
        'observatory.LS240_ID.feeds.temperatures.Channel 7 T',
        'observatory.LS240_ID.feeds.temperatures.Channel 5 T',
    ]

* **Option 3**

Use a configuration file. See Below.

Define How the data is returned
``````````````````````````````````````````````

The data is returned as a dictionary of the format::

    {
        'name' : (time, data)
    }

``time`` and ``data`` are arrays of the times / data from each loaded field

* **Option 1** 

No keyword arguments means ``load_range`` will return ``name`` set to be the
field name. But this is long.

* **Option 2** 

Use the ``alias = [list, of, desired, names]`` which must be the length of 
``fields``. Now the dictionary will have these alias as the ``name``.

* **Option 3**

Use a configuration file. See Below.

Create a Configuration file
``````````````````````````````````````````````

Because why deal with all these individual pieces when you don't have to?

Define a ``yaml`` file and pass it to ``load_range`` with the ``config``
keyword argument. Ex File::

    data_dir: '/data/ocs'
    
    field_list:
        '40k_dr_side' : 'observatory.LEIA.feeds.temperatures.Channel 7 T'
        '40k_far_side': 'observatory.LEIA.feeds.temperatures.Channel 8 T'
        '80k_dr_side' : 'observatory.LEIA.feeds.temperatures.Channel 5 T' 
        '80k_far_side': 'observatory.LEIA.feeds.temperatures.Channel 6 T'
        '4k_far_side' : 'observatory.YODA.feeds.temperatures.Channel 1 T'
        '4k_dr_side'  : 'observatory.YODA.feeds.temperatures.Channel 2 T'

``data_dir`` sets the directory and ``field_list`` has the list of ``'alias':'field'``.


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


Low-level HK Stream Processing
==============================

There are a few HK stream processing objects that are intended for use
as modules in a ``G3Pipeline``.  These are their stories.

Module: HKScanner
-----------------

.. autoclass:: so3g.hk.HKScanner
   :members:

Module: HKReframer
------------------

.. autoclass:: so3g.hk.HKReframer
   :members:


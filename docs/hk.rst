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
-------------------------------

One of the most basic things we might want to do is load data between
a time range. For `.g3` files that are saved by an OCS Aggregator, there
is a specific folder structure and file naming scheme. :func:`so3g.hk.load_range`
is writen to load data for a specified time frame saved in that style.

Example Use::

    from so3g.hk import load_range    
    data = load_range(start, stop, **kwargs)

Defining Time Ranges
````````````````````

There are several options for defining start and stop. These options
are passed to :func:`so3g.hk.getdata.to_timestamp` to be parsed.

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
`````````````````````````````

* **Option 1** 

  Use the ``data_dir='/path/to/ocs/hk/data'`` keyword argument. This
  should be the directory where all the first five digit ctime folders
  are located. ``load_range`` will look for the data there.

* **Option 2** 

  Set an environment variable ``export
  OCS_DATA_DIR=/path/to/ocs/hk/data``.  ``load_range`` will will
  automatically look there if it isn't overridden by option 1.

* **Option 3**

  Use a configuration file. See Below.


Define Which Data to Load
`````````````````````````

* **Option 1** 

  No keyword arguments means ``load_range`` will return *every* field
  it can find. this might take a very long time.

* **Option 2** 

  Use the ``fields = [list, of, field, names]`` keyword argument. Example::
    
    fields = [
        'observatory.LS240_ID.feeds.temperatures.Channel_7_T',
        'observatory.LS240_ID.feeds.temperatures.Channel_5_T',
    ]

* **Option 3**

  Use a configuration file. See Below.


Define How the Data is Returned
```````````````````````````````

The data is returned as a dictionary of the format::

    {
        'name' : (time, data)
    }

``time`` and ``data`` are arrays of the times / data from each loaded field

* **Option 1** 

  No keyword arguments means ``load_range`` will return ``name`` set
  to be the field name. But this is long.

* **Option 2** 

  Use the ``alias = [list, of, desired, names]`` which must be the
  length of ``fields``. Now the dictionary will have these alias as
  the ``name``.

* **Option 3**

  Use a configuration file. See Below.


Create a Configuration file
```````````````````````````

Because why deal with all these individual pieces when you don't have to?

Create a ``yaml`` file and pass the filename to ``load_range`` with
the ``config`` keyword argument. Example file::

    data_dir: '/data/ocs'
    
    field_list:
        '40k_dr_side' : 'observatory.LEIA.feeds.temperatures.Channel_7_T'
        '40k_far_side': 'observatory.LEIA.feeds.temperatures.Channel_8_T'
        '80k_dr_side' : 'observatory.LEIA.feeds.temperatures.Channel_5_T'
        '80k_far_side': 'observatory.LEIA.feeds.temperatures.Channel_6_T'
        '4k_far_side' : 'observatory.YODA.feeds.temperatures.Channel_1_T'
        '4k_dr_side'  : 'observatory.YODA.feeds.temperatures.Channel_2_T'

``data_dir`` sets the directory and ``field_list`` has the list of ``'alias':'field'``.


Function References
```````````````````

.. autofunction:: so3g.hk.load_range

.. autofunction:: so3g.hk.getdata.to_timestamp


Exploring HK Data
-----------------

The HKTree object provides a way to browse through available HK data
fields from within an interactive python session (such as through the
python or ipython interpreters or in a jupyter session).

Instead of accessing HK fields through their long string names, such as::

  field_name = "observatory.hk.rotation.feeds.hwprotation.kikusui_curr"

a field can be referred to as a named object in a hierarchy of
attributes::

  tree = HKTree('2022-10-20', '2022-10-22', config='hk.yaml')
  field = tree.hk.rotation.hwprotation.kikusui_curr

By exposing the available fields as a tree of attributes, rather than
huge set of long string keys, a user working interactively can use
tab-completion to find fields easily.

Having identified a field or set of fields of interest, the data can
be loaded by calling pseudo-private methods on the field directly, e.g.::

  data_dict = field._load()

After calling load, the data are stored in the field reference itself,
and so can be retrieved with::

  times, values = field._data

The ``_load`` method can be called on "non-terminal" nodes of the
tree, for example::

  data_dict = tree.hk.rotation._load()

would load all the fields that are sub- (or sub-sub-, ...) attributes
of ``tree.hk.rotation``.

The system is intended to complement :func:`so3g.hk.load_range`, and
uses the same sort of configuration file.

See more detailed examples below, as well as the :ref:`Class Reference
<HKTree Class Reference>`.

Instantiating an HKTree
```````````````````````

To create an HKTree requires at least a path to an "OCS style" data
directory (in the same sense as :func:`so3g.hk.load_range`)::

  tree = hk.HKTree(data_dir='/mnt/so1/data/ucsd-sat1/hk/')

By default, the returned object will only look through data from the
past 24 hours.  To specify a range of dates of interest, use the start
and stop parameters::

  tree = hk.HKTree(data_dir='/mnt/so1/data/ucsd-sat1/hk/',
                   start='2022-10-20', stop='2022-10-22')

As with ``load_range``, passing ``data_dir`` is not necessary if the
``OCS_DATA_DIR`` environment variable is set.  However, a config file
may be the best way to go.


Configuration file
``````````````````

The configuration file syntax is as in ``load_range``.  However you
can also specify:

- ``pre_proc_dir``: The default value for ``pre_proc_dir``.
- ``skip_tokens``: Tokens to append to the ``skip`` parameter for
  ``HKTree()``.

Aliases are treated in a special way by HKTree; see `Using ._aliases`
below.


Finding and Loading data
````````````````````````

The nodes in the tree are all either "terminal" or not.  The terminal
nodes are associated with a single specific HK field
(e.g. ``tree.hk.rotation.hwprotation.kikusui_curr``) while
non-terminal nodes are do not have associated fields, but have child
attributes (e.g. ``tree.hk.rotation``).

Data is loaded by calling ``._load()`` on any node in the tree.  This
function will return a dictionary of data, mapping the full field names
to tuples of data (this is similar to what ``load_range`` returns).

Following load, the tuples are also available in the ``._data``
attribute of each terminal node.  For example::

  data_dict = tree.hk.rotation._load()

will return a dict with multiple fields.  But after that call one can
also access single field data on terminal nodes, e.g.::

  t, val = tree.hk.rotation.kikusui_curr

If you want to clear RAM / start over, call ``._clear()`` on any node
in the tree to clear the data from all its child nodes.  E.g.::

  tree.hk._clear()


Using ._aliases
```````````````

The ``load_range`` function permits users to associate aliases with
long field names, using the ``field_list`` in the config file (or the
``aliases`` argument).  The ``HKTree`` makes those fields available,
under their aliases, in a special attribute called aliases.  For
example the config assignments like this::

    field_list:
        't40k_dr_side' : 'observatory.LEIA.feeds.temperatures.Channel_7_T'
        't40k_far_side': 'observatory.LEIA.feeds.temperatures.Channel_8_T'

would lead to the existence of attributes::

  tree._aliases.t40k_dr_side
  tree._aliases.t40k_far_side

(Note that attributes can't be accessed if they begin with a digit or
contain special characters... so avoid them in your field_list.)

Fields exposed under ``._aliases`` also exist in the full tree -- the
terminal attributes here are actually the same ones .  You
can run ``._load()`` and ``._clear()`` on ``._aliases`` and it will
operate on all the attributes contained there.

In particular note that ``tree._aliases._load()`` should return a data
dictionary where the alias names are used as the keys, instead of the

You can dynamicaly add new fields to the aliases list (as a way of
grouping things together under shorter names) using
``tree._add_alias()``, providing an alias string and a target field.
The following are equivalent::

  tree._add_alias('short_name', tree.LEIA.temperatures.Channel_1_T)
  tree._add_alias('short_name', 'observatory.tree.LEIA.feeds.temperatures.Channel_1_T')

You can add a sub-tree with children; the aliases will be generated
automatically; for example:

  tree._add_alias('terms', tree.LEIA.temperatures)

Would create aliases called ``'therms_Channel_1T'``,
``'therms_Channel_2T'``, etc.


.. _HKTree Class Reference:

Class Reference
```````````````

You should see auto-generated class documentation below.

.. autoclass:: so3g.hk.tree.HKTree
   :members: _load, _clear, _add_alias

.. autoclass:: so3g.hk.tree.HKRef
   :members: _load, _clear


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

Checking Files with so-hk-tool
------------------------------

The command line tool ``so-hk-tool`` can be used to scan one or more SO
HK files and summarize info about the files, the providers, and the
fields within.  For example::

  $ so-hk-tool list-provs /mnt/so1/data/ucsd-sat1/hk/16685/1668577223.g3
  provider_name                           total_bytes frame_bytes
  --------------------------------------- ----------- -----------
  BK9130C-1.psu_output                          67152     11192.0
  BK9130C-2.psu_output                          67152     11192.0
  DRptc1.ptc_status                           1006500     16775.0
  LSA21US.temperatures                         388271      6365.1
  LSA21YC.temperatures                         869316    144886.0
  LSA22QC.temperatures                         585204     97534.0
  LSA24LY.temperatures                         897660     14961.0
  LSA24M5.temperatures                         896376     14939.6
  LSA2619.temperatures                         895620     14927.0
  SSG3KRM-2_2.ups                               34169      5694.8
  ...

The tool can run the following analyses:

- ``list-files``: List each file, its size, and report any stray bytes
  (partial trailing frames).
- ``list-provs``: List each provider found in the files.
- ``list-fields``: List each field found in the files.

See more details below.  Please note that when presenting provider and
field names, the program strips out tokens ``observatory`` and
``feeds``, by default (for example, ``observatory.DRptc1.feeds.ptc_status``
becomes ``DRptc1.ptc_status``).  Pass ``--strip-tokens=""`` to instead
show the full provider / feed names.

.. argparse::
   :prog: so-hk-tool
   :module: so3g.hk.cli
   :func: get_parser


HK Data Types and File Structure
--------------------------------

The HK file structures and versions are described in
https://github.com/simonsobs/tod2maps_docs/tod_format.

As of August 2021, all HK data uses schema version 2, which supports
vectors of float, integer, or string data.  The original form (schema
version 0) supported only floats, and required access to compiled G3
extensions in so3g.  In case there are still some schema 0 data out
there, conversion of individual files can be achieved with help from
:class:`so3g.hk.HKTranslator`.


Writing HK Data
---------------

The so3g.hk module provides limited assistance with creating HK data
files.  The :class:`so3g.hk.HKSessionHelper` may be used to produce
template frames that can be used as a basis for an HK data stream.
However, the code in this module does not enforce validity.  (The OCS
"aggregator" Agent has more sophisticated logic to help write only
valid HK frame streams.)

Here is a short example that creates a housekeeping file containing
some fake pointing information (it can be found the repository as
``demos/write_hk.py``):

.. include:: ../demos/write_hk.py
    :literal:

When extending this example for other purposes, here are a few things
to remember, to help generate *valid* HK streams:

- Notice that "block" is a G3TimesampleMap, a class designed to store
  multiple data vectors alongside a single vector of timesamples.  If
  your provider has fields with different time sampling, group them so
  each block corresponds to mutually co-sampled fields.
- The "block name" is an internal bookkeeping thing and won't be
  visible to the consumer of the HK data.  In the example above, the
  data vectors would be exposed through the combination of the
  provider and field name, i.e.. "platform.az", "platform.el".
- Once you start a provider, and get a ``prov_id`` for it, then any
  named blocks you add must always have the same fields with the same
  data types.  (For example, it would be illegal to only have some
  frames where the block contains only "el" and not "az".)  If you
  need to change the list of fields, use ``remove_provider()`` and the
  re-add the provider (with the same name).

Below, find the documentation for the HKSessionHelper class.

.. autoclass:: so3g.hk.HKSessionHelper
   :members:


Low-level HK Stream Processing
------------------------------

There are a few HK stream processing objects that are intended for use
as modules in a ``G3Pipeline``.  These are their stories.

Module: HKScanner
`````````````````

.. autoclass:: so3g.hk.HKScanner
   :members:

Module: HKReframer
``````````````````

.. autoclass:: so3g.hk.HKReframer
   :members:

Module: HKTranslator
````````````````````

.. autoclass:: so3g.hk.HKTranslator
   :members:


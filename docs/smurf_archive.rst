Loading Detector Data
======================

The SmurfArchive module contains code for indexing and quickly loading detector
data from an archive of streamed data written by the smurf-recorder OCS agent.
A `data archive` should contain all timestream data recorded on a given system.
The file paths within a data archive will look like::

    <archive-path>/<date-code>/<stream-id>/<filename>.g3

For each stream-id, it is important that the time ranges of different files do
not overlap or else data may not be returned in the correct order. This current
version of the SmurfArchive system only works for a single stream-id, however
future changes to the SmurfRecorder OCS agent and the SmurfArchive class will
provide support for archives containing timestreams from multiple stream-id's.
This change will be necessary for systems with multiple smurf cards.


Indexing Timestream Data
--------------------------
To index data, the SmurfArchive module will walk through an archive and record
information about each frame to an sqlite database, including the frame time,
frame type, the number of channels and samples, etc. Only one index database
is needed for each archive. The following snippet will write an index database
to ``<archive_path>/frames.db`` or update an existing database with new frames
if it already exists::

    from so3g.smurf import SmurfArchive
    arc = SmurfArchive('/data/timestreams')
    arc.index_archive()

In this example, all detector data stored in ``/data/timestreams`` will be
indexed and the index database will be written to
``/data/timestreams/frames.db``.

A different database path can be specified by passing the ``db_path`` argument
to the ``SmurfArchive``.

Loading Data
-------------

To load detector data that includes the range ``[start, stop]``, use::

    from so3g.smurf import SmurfArchive
    arc = SmurfArchive('/data/timestreams')
    d = arc.load_data(start, end)

The data will be loaded into a named_tuple which contains the detector data,
primary data, and TES bias values, and metadata status. See the `load_data` Api
for a more detailed description of allowed keyword arguments and the return
value.

API
---

.. autoclass:: so3g.smurf.SmurfArchive.SmurfArchive
   :members: index_archive, load_data, load_status


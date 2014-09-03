.. _creating_frontend:

Creating A New Code Frontend
============================

.. warning: This section is not yet updated to work with yt 3.0.  If you
            have a question about making a custom derived quantity, please
            contact the mailing list.

yt is designed to support analysis and visualization of data from multiple
different simulation codes, although it has so far been most successfully
applied to Adaptive Mesh Refinement (AMR) data. For a list of codes and the
level of support they enjoy, see :ref:`code-support`.

We'd like to support a broad range of codes, both AMR-based and otherwise. To
add support for a new code, a few things need to be put into place. These
necessary structures can be classified into a couple categories:

 * Data meaning: This is the set of parameters that convert the data into
   physically relevant units; things like spatial and mass conversions, time
   units, and so on.
 * Data localization: These are structures that help make a "first pass" at data
   loading. Essentially, we need to be able to make a first pass at guessing
   where data in a given physical region would be located on disk. With AMR
   data, this is typically quite easy: the grid patches are the "first pass" at
   localization.
 * Data reading: This is the set of routines that actually perform a read of
   either all data in a region or a subset of that data.

Data Meaning Structures
-----------------------

If you are interested in adding a new code, be sure to drop us a line on
`yt-dev <http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org>`_!

To get started, make a new directory in ``yt/frontends`` with the name of your
code -- you can start by copying into it the contents of the ``stream``
directory, which is a pretty empty format. You'll then have to create a subclass
of ``Dataset``. This subclass will need to handle conversion between the
different physical units and the code units; for the most part, the examples of
``OrionDataset`` and ``EnzoDataset`` should be followed, but
``ChomboDataset``, as a slightly newer addition, can also be used as an
instructive example -- be sure to add an ``_is_valid`` classmethod that will
verify if a filename is valid for that output type, as that is how "load" works.

A new set of fields must be added in the file ``fields.py`` in that directory.
For the most part this means subclassing ``CodeFieldInfoContainer`` and adding
the necessary fields specific to that code. Here is the Chombo field container:

.. code-block:: python

    from UniversalFields import *
    class ChomboFieldContainer(CodeFieldInfoContainer):
        _shared_state = {}
        _field_list = {}
    ChomboFieldInfo = ChomboFieldContainer()
    add_chombo_field = ChomboFieldInfo.add_field

The field container is a shared state object, which is why we explicitly set
``_shared_state`` equal to a mutable.

Data Localization Structures
----------------------------

As of right now, the "grid patch" mechanism is going to remain in yt, however in
the future that may change. As such, some other output formats -- like Gadget --
may be shoe-horned in, slightly.

Hierarchy
^^^^^^^^^

To set up data localization, an ``AMRHierarchy`` subclass must be added in the
file ``data_structures.py``. The index object must override the following
methods:

 * ``_detect_fields``: ``self.field_list`` must be populated as a list of
   strings corresponding to "native" fields in the data files.
 * ``_setup_classes``: it's probably safe to crib this from one of the other
   ``AMRHierarchy`` subclasses.
 * ``_count_grids``: this must set self.num_grids to be the total number of
   grids in the simulation.
 * ``_parse_index``: this must fill in ``grid_left_edge``,
   ``grid_right_edge``, ``grid_particle_count``, ``grid_dimensions`` and
   ``grid_levels`` with the appropriate information. Additionally, ``grids``
   must be an array of grid objects that already know their IDs.
 * ``_populate_grid_objects``: this initializes the grids by calling
   ``_prepare_grid`` and ``_setup_dx`` on all of them.  Additionally, it should
   set up ``Children`` and ``Parent`` lists on each grid object.
 * ``_setup_unknown_fields``: If a field is in the data file that yt doesn't
   already know, this is where you make a guess at it.
 * ``_setup_derived_fields``: ``self.derived_field_list`` needs to be made a
   list of strings that correspond to all derived fields valid for this
   index.

For the most part, the ``ChomboHierarchy`` should be the first place to look for
hints on how to do this; ``EnzoHierarchy`` is also instructive.

Grids
^^^^^

A new grid object, subclassing ``AMRGridPatch``, will also have to be added.
This should go in ``data_structures.py``. For the most part, this may be all
that is needed:

.. code-block:: python

    class ChomboGrid(AMRGridPatch):
        _id_offset = 0
        __slots__ = ["_level_id"]
        def __init__(self, id, index, level = -1):
            AMRGridPatch.__init__(self, id, filename = index.index_filename,
                                  index = index)
            self.Parent = []
            self.Children = []
            self.Level = level


Even the most complex grid object, ``OrionGrid``, is still relatively simple.

Data Reading Functions
----------------------

In ``io.py``, there are a number of IO handlers that handle the mechanisms by
which data is read off disk.  To implement a new data reader, you must subclass
``BaseIOHandler`` and override the following methods:

 * ``_read_field_names``: this routine accepts a grid object and must return all
   the fields in the data file affiliated with that grid. It is used at the
   initialization of the ``AMRHierarchy`` but likely not later.
 * ``modify``: This accepts a field from a data file and returns it ready to be
   used by yt. This is used in Enzo data for preloading.
 * ``_read_data_set``: This accepts a grid object and a field name and must
   return that field, ready to be used by yt as a NumPy array. Note that this
   presupposes that any actions done in ``modify`` (above) have been executed.
 * ``_read_data_slice``: This accepts a grid object, a field name, an axis and
   an (integer) coordinate, and it must return a slice through the array at that
   value.
 * ``preload``: (optional) This accepts a list of grids and a list of datasets
   and it populates ``self.queue`` (a dict keyed by grid id) with dicts of
   datasets.
 * ``_read_exception``: (property) This is a tuple of exceptions that can be
   raised by the data reading to indicate a field does not exist in the file.


And that just about covers it. Please feel free to email
`yt-users <http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org>`_ or
`yt-dev <http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org>`_ with
any questions, or to let us know you're thinking about adding a new code to yt.

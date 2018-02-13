.. _creating_frontend:

Creating A New Code Frontend
============================

.. warning: This section is not yet updated to work with yt 3.0.  If you
            have a question about making a custom derived quantity, please
            contact the mailing list.

yt is designed to support analysis and visualization of data from
multiple different simulation codes. For a list of codes and the level
of support they enjoy, see :ref:`code-support`.

We'd like to support a broad range of codes, both Adaptive Mesh
Refinement (AMR)-based and otherwise. To add support for a new code, a
few things need to be put into place. These necessary structures can
be classified into a couple categories:

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
`yt-dev <https://mail.python.org/mm3/archives/list/yt-dev@python.org/>`_!

To get started, make a new directory in ``yt/frontends`` with the name
of your code and add the name into ``yt/frontends/api.py``.
Copying the contents of the ``yt/frontends/_skeleton``
directory will add a lot of boilerplate for the required classes and
methods that are needed.  In particular, you'll have to create a
subclass of ``Dataset`` in the data_structures.py file. This subclass
will need to handle conversion between the different physical units
and the code units (typically in the ``_set_code_unit_attributes()``
method), read in metadata describing the overall data on disk (via the
``_parse_parameter_file()`` method), and provide a ``classmethod``
called ``_is_valid()`` that lets the ``yt.load`` method help identify an
input file as belonging to *this* particular ``Dataset`` subclass.
For the most part, the examples of
``yt.frontends.boxlib.data_structures.OrionDataset`` and
``yt.frontends.enzo.data_structures.EnzoDataset`` should be followed,
but ``yt.frontends.chombo.data_structures.ChomboDataset``, as a
slightly newer addition, can also be used as an instructive example.

A new set of fields must be added in the file ``fields.py`` in your
new directory.  For the most part this means subclassing
``FieldInfoContainer`` and adding the necessary fields specific to
your code. Here is a snippet from the base BoxLib field container:

.. code-block:: python

    from yt.fields.field_info_container import FieldInfoContainer
    class BoxlibFieldInfo(FieldInfoContainer):
        known_other_fields = (
            ("density", (rho_units, ["density"], None)),
	    ("eden", (eden_units, ["energy_density"], None)),
	    ("xmom", (mom_units, ["momentum_x"], None)),
	    ("ymom", (mom_units, ["momentum_y"], None)),
	    ("zmom", (mom_units, ["momentum_z"], None)),
	    ("temperature", ("K", ["temperature"], None)),
	    ("Temp", ("K", ["temperature"], None)),
	    ("x_velocity", ("cm/s", ["velocity_x"], None)),
	    ("y_velocity", ("cm/s", ["velocity_y"], None)),
	    ("z_velocity", ("cm/s", ["velocity_z"], None)),
	    ("xvel", ("cm/s", ["velocity_x"], None)),
	    ("yvel", ("cm/s", ["velocity_y"], None)),
	    ("zvel", ("cm/s", ["velocity_z"], None)),
	)

	known_particle_fields = (
	    ("particle_mass", ("code_mass", [], None)),
	    ("particle_position_x", ("code_length", [], None)),
	    ("particle_position_y", ("code_length", [], None)),
	    ("particle_position_z", ("code_length", [], None)),
	    ("particle_momentum_x", (mom_units, [], None)),
	    ("particle_momentum_y", (mom_units, [], None)),
	    ("particle_momentum_z", (mom_units, [], None)),
	    ("particle_angmomen_x", ("code_length**2/code_time", [], None)),
	    ("particle_angmomen_y", ("code_length**2/code_time", [], None)),
	    ("particle_angmomen_z", ("code_length**2/code_time", [], None)),
	    ("particle_id", ("", ["particle_index"], None)),
	    ("particle_mdot", ("code_mass/code_time", [], None)),
	)

The tuples, ``known_other_fields`` and ``known_particle_fields``
contain entries, which are tuples of the form ``("name", ("units",
["fields", "to", "alias"], "display_name"))``.  ``"name"`` is the name
of a field stored on-disk in the dataset. ``"units"`` corresponds to
the units of that field.  The list ``["fields", "to", "alias"]``
allows you to specify additional aliases to this particular field; for
example, if your on-disk field for the x-direction velocity were
``"x-direction-velocity"``, maybe you'd prefer to alias to the more
terse name of ``"xvel"``.  ``"display_name"`` is an optional parameter
that can be used to specify how you want the field to be displayed on
a plot; this can be LaTeX code, for example the density field could
have a display name of ``r"\rho"``.  Omitting the ``"display_name"``
will result in using a capitalized version of the ``"name"``.

.. _bfields-frontend:

Creating Aliases for Magnetic Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setting up access to the magnetic fields in your dataset requires special
handling, because in different unit systems magnetic fields have different
dimensions (see :ref:`bfields` for an explanation). If your dataset includes
magnetic fields, you should include them in ``known_other_fields``, but do
not set up aliases for them--instead use the special handling function
:meth:`~yt.fields.magnetic_fields.setup_magnetic_field_aliases`. It takes
as arguments the ``FieldInfoContainer`` instance, the field type of the
frontend, and the list of magnetic fields from the frontend. Here is an
example of how this is implemented in the FLASH frontend:

.. code-block:: python

    class FLASHFieldInfo(FieldInfoContainer):
        known_other_fields = (
            ...
            ("magx", (b_units, [], "B_x")), # Note there is no alias here
            ("magy", (b_units, [], "B_y")),
            ("magz", (b_units, [], "B_z")),
            ...
        )

        def setup_fluid_fields(self):
            from yt.fields.magnetic_field import \
                setup_magnetic_field_aliases
            ...
            setup_magnetic_field_aliases(self, "flash", ["mag%s" % ax for ax in "xyz"])

This function should always be imported and called from within the
``setup_fluid_fields`` method of the ``FieldInfoContainer``. If this
function is used, converting between magnetic fields in different
:ref:`unit systems <unit_systems>` will be handled automatically.

Data Localization Structures
----------------------------

These functions and classes let yt know about how the arrangement of
data on disk corresponds to the physical arrangement of data within
the simulation.  yt has grid datastructures for handling both
patch-based and octree-based AMR codes.  The terms 'patch-based'
and 'octree-based' are used somewhat loosely here.  For example,
traditionally, the FLASH code used the paramesh AMR library, which is
based on a tree structure, but the FLASH frontend in yt utilizes yt's
patch-based datastructures.  It is up to the frontend developer to
determine which yt datastructures best match the datastructures of
their simulation code.

Both approaches -- patch-based and octree-based -- have a concept of a
*Hierarchy* or *Index* (used somewhat interchangeably in the code) of
datastructures and something that describes the elements that make up
the Hierarchy or Index.  For patch-based codes, the Index is a
collection of ``AMRGridPatch`` objects that describe a block of zones.
For octree-based codes, the Index contains datastructures that hold
information about the individual octs, namely an ``OctreeContainer``.

Hierarchy or Index
^^^^^^^^^^^^^^^^^^

To set up data localization, a ``GridIndex`` subclass for patch-based
codes or an ``OctreeIndex`` subclass for octree-based codes must be
added in the file ``data_structures.py``. Examples of these different
types of ``Index`` can be found in, for example, the
``yt.frontends.chombo.data_structures.ChomboHierarchy`` for patch-based
codes and ``yt.frontends.ramses.data_structures.RAMSESIndex`` for
octree-based codes.

For the most part, the ``GridIndex`` subclass must override (at a
minimum) the following methods:

 * ``_detect_output_fields()``: ``self.field_list`` must be populated as a list
   of strings corresponding to "native" fields in the data files.
 * ``_count_grids()``: this must set ``self.num_grids`` to be the total number
   of grids (equivalently ``AMRGridPatch``'es) in the simulation.
 * ``_parse_index()``: this must fill in ``grid_left_edge``,
   ``grid_right_edge``, ``grid_particle_count``, ``grid_dimensions`` and
   ``grid_levels`` with the appropriate information.  Each of these variables
   is an array, with an entry for each of the ``self.num_grids`` grids.
   Additionally, ``grids``  must be an array of ``AMRGridPatch`` objects that
   already know their IDs.
 * ``_populate_grid_objects()``: this initializes the grids by calling
   ``_prepare_grid()`` and ``_setup_dx()`` on all of them.  Additionally, it
   should set up ``Children`` and ``Parent`` lists on each grid object.

The ``OctreeIndex`` has somewhat analogous methods, but often with
different names; both ``OctreeIndex`` and ``GridIndex`` are subclasses
of the ``Index`` class.  In particular, for the ``OctreeIndex``, the
method ``_initialize_oct_handler()`` setups up much of the oct
metadata that is analogous to the grid metadata created in the
``GridIndex`` methods ``_count_grids()``, ``_parse_index()``, and
``_populate_grid_objects()``.

Grids
^^^^^

.. note:: This section only applies to the approach using yt's patch-based
	  datastructures.  For the octree-based approach, one does not create
	  a grid object, but rather an ``OctreeSubset``, which has methods
	  for filling out portions of the octree structure.  Again, see the
	  code in ``yt.frontends.ramses.data_structures`` for an example of
	  the octree approach.

A new grid object, subclassing ``AMRGridPatch``, will also have to be added in
``data_structures.py``. For the most part, this may be all
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


Even one of the more complex grid objects,
``yt.frontends.boxlib.BoxlibGrid``, is still relatively simple.

Data Reading Functions
----------------------

In ``io.py``, there are a number of IO handlers that handle the
mechanisms by which data is read off disk.  To implement a new data
reader, you must subclass ``BaseIOHandler``.  The various frontend IO
handlers are stored in an IO registry - essentially a dictionary that
uses the name of the frontend as a key, and the specific IO handler as
a value.  It is important, therefore, to set the ``dataset_type``
attribute of your subclass, which is what is used as the key in the IO
registry.  For example:

.. code-block:: python

    class IOHandlerBoxlib(BaseIOHandler):
        _dataset_type = "boxlib_native"
	...

At a minimum, one should also override the following methods

* ``_read_fluid_selection()``: this receives a collection of data "chunks", a
  selector describing which "chunks" you are concerned with, a list of fields,
  and the size of the data to read.  It should create and return a dictionary
  whose keys are the fields, and whose values are numpy arrays containing the
  data.  The data should actually be read via the ``_read_chunk_data()``
  method.
* ``_read_chunk_data()``: this method receives a "chunk" of data along with a
  list of fields we want to read.  It loops over all the grid objects within
  the "chunk" of data and reads from disk the specific fields, returning a
  dictionary whose keys are the fields and whose values are numpy arrays of
  the data.

If your dataset has particle information, you'll want to override the
``_read_particle_coords()`` and ``read_particle_fields()`` methods as
well.  Each code is going to read data from disk in a different
fashion, but the ``yt.frontends.boxlib.io.IOHandlerBoxlib`` is a
decent place to start.

And that just about covers it. Please feel free to email
`yt-users <https://mail.python.org/mm3/archives/list/yt-users@python.org/>`_ or
`yt-dev <https://mail.python.org/mm3/archives/list/yt-dev@python.org/>`_ with
any questions, or to let us know you're thinking about adding a new code to yt.

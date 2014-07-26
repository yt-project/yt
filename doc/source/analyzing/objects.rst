XXX.. _using-objects:

Data Objects
============

What are Data Objects in yt?
----------------------------

Data objects (also called *Data Containers*) are used in yt as convenience 
structures for grouping data in logical ways that make sense in the context 
of the dataset as a whole.  Some of the data objects are geometrical groupings 
of data (e.g. sphere, region--a 3D box, cylinder, etc.).  Others represent 
data products derived from your dataset (e.g. slices, streamlines, surfaces).
Still other data objects group multiple objects together or filter them
(e.g. data dollection, cut region).  

To generate standard plots, objects rarely need to be directly constructed.
However, for detailed data inspection as well as hand-crafted derived data,
objects can be exceptionally useful and even necessary.

For geometric objects, if the shape intersectsXXX

How to Create an Object
-----------------------

To create an object, you usually only need a loaded dataset, the name of 
the object type, and the relevant parameters for your object.  Here is a common
example for creating a ``Region`` object that covers all of your data volume.

.. code-block:: python

   import yt
   ds = yt.load("RedshiftOutput0005")
   ad = ds.all_data()

Alternatively, we could create a sphere object of radius 1 kpc on location 
[0.5, 0.5, 0.5] using the dataset quantity 1 kpc:

.. code-block:: python

   import yt
   ds = yt.load("RedshiftOutput0005")
   sp = ds.sphere([0.5, 0.5, 0.5], ds.quan(1, 'kpc'))

.. _available-objects:

Available Objects
-----------------

As noted above, there are numerous types of objects.  Here we group them
into:

* *Geometric Objects* - Data is selected based on spatial shapes in the dataset
* *Filtering Objects* - Data is selected based on other field criteria
* *Collection Objects* - Multiple objects grouped together
* 

Geometric Objects
^^^^^^^^^^^^^^^^^

0D
""

**Point** 
    Aliased to :class:`~yt.data_objects.data_containers.YTPointBase`    
    Usage: ``point(coords)``
    A zero-dimensional point defined by a single cell at specified coordinates.

1D
""

**Axis-Aligned Ray** (aliased to :class:`~yt.data_objects.data_containers.YTOrthoRayBase`)
    | Usage: ``ortho_ray()``
    | A one-dimensional line of data cells stretching through the full domain aligned with one of the x,y,z axes.

**Arbitrary-Aligned Ray** (aliased to :class:`~yt.data_objects.data_containers.YTRayBase`)
    | Usage: ``ray()``
    | A one-dimensional line of data cells stretching through the full domain defined by arbitrary start and end coordinates.

2D 
""

**Axis-Aligned Slice** (aliased to :class:`~yt.data_objects.data_containers.YTSliceBase`)
    | Usage: ``slice()``

**Arbitrary-Aligned Slice** (aliased to :class:`~yt.data_objects.data_containers.YTCuttingPlaneBase`)
    | Usage: ``cutting()``

3D
""

**Disk/Cylinder** (aliased to :class:`~yt.data_objects.data_containers.YTDiskBase`)
    | Usage: ``disk()``

**Box Region** (aliased to :class:`~yt.data_objects.data_containers.YTRegionBase`)
    | Usage: ``region()``

**Sphere** (aliased to :class:`~yt.data_objects.data_containers.YTSphereBase`)
    | Usage: ``sphere()``

**Ellipsoid** (aliased to :class:`~yt.data_objects.data_containers.YTEllipsoidBase`)
    | Usage: ``ellipsoid()``

**All Data** (aliased to :class:`~yt.data_objects.data_containers.YTRegionBase`)
    | Usage: ``all_data()``

Filtering and Grouping Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Boolean Regions** (Note: not yet implemented in yt 3.0)
    | Usage: ``boolean()``

**Mesh Field Filter** (aliased to :class:`~yt.data_objects.data_containers.YTCutRegionBase`)
    | Usage: ``cut_region()``

**Collection of Data Objects** (aliased to :class:`~yt.data_objects.data_containers.YTDataCollectionBase`)
    | Usage: ``data_collection()``

Data Product Objects
^^^^^^^^^^^^^^^^^^^^

**Streamline** (aliased to :class:`~yt.data_objects.data_containers.YTStreamlineBase`)
    | Usage: ``streamline()``

**Projection** (aliased to :class:`~yt.data_objects.data_containers.YTQuadTreeProjBase`)
    | Usage: ``proj()``

**Fixed-Resolution Region** (aliased to :class:`~yt.data_objects.data_containers.YTCoveringGridBase`)
    | Usage: ``covering_grid()``

**Fixed-Resolution Region with Smoothing** (aliased to :class:`~yt.data_objects.data_containers.YTSmoothedCoveringGridBase`)
    | Usage: ``smoothed_covering_grid()``

**Fixed-Resolution Region for Particle Deposition** (aliased to :class:`~yt.data_objects.data_containers.YTArbitraryGridBase `)
    | Usage: ``arbitrary_grid()``

**Surface** (aliased to :class:`~yt.data_objects.data_containers.YTSurfaceBase`)
    | Usage: ``surface()``



The following objects are available, all of which hang off of the index
object.  To access them, you would do something like this (as for a
:class:`region`):

.. code-block:: python

   import yt
   ds = yt.load("RedshiftOutput0005")
   reg = ds.region([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])


.. _boolean_data_objects:

Combining Objects: Boolean Data Objects
---------------------------------------

.. note:: Boolean Data Objects have not yet been ported to yt 3.0 from
    yt 2.x.  If you are interested in aiding in this port, please contact
    the yt-dev mailing list.  Until it is ported, this functionality below
    will not work.

A special type of data object is the *boolean* data object.
It works only on three-dimensional objects.
It is built by relating already existing data objects with boolean operators.
The boolean logic may be nested using parentheses, and
it supports the standard "AND", "OR", and "NOT" operators:

* **"AND"** When two data objects are related with an "AND", the combined
  data object is the volume of the simulation covered by both objects, and
  not by just a single object.
* **"OR"** When two data objects are related with an "OR", the combined
  data object is the volume(s) of the simulation covered by either of the
  objects.
  For example, this may be used to combine disjoint objects into one.
* **"NOT"** When two data objects are related with a "NOT", the combined
  data object is the volume of the first object that the second does not
  cover.
  For example, this may be used to cut out part(s) of the first data object
  utilizing the second data object.
* **"(" or ")"** Nested logic is surrounded by parentheses. The order of
  operations is such that the boolean logic is evaluated inside the
  inner-most parentheses, first, then goes upwards.
  The logic is read left-to-right at all levels (crucial for the "NOT"
  operator).

Please see the :ref:`cookbook` for some examples of how to use the boolean
data object.

.. _derived-quantities:

Processing Objects: Derived Quantities
--------------------------------------

Derived quantities are a way of operating on a collection of cells and
returning a set of values that is fewer in number than the number of cells --
yt already knows about several.  Every 3D data object (see
:ref:`using-objects`) provides a mechanism for access to derived quantities.
These can be accessed via the ``quantities`` interface, like so:

.. code-block:: python

   ds = load("my_data")
   dd = ds.all_data()
   dd.quantities.angular_momentum_vector()

The following quantities are available via the ``quantities`` interface.

.. include:: _dq_docstrings.inc

Creating Derived Quantities
+++++++++++++++++++++++++++

The basic idea is that you need to be able to operate both on a set of data,
and a set of sets of data.  (If this is not possible, the quantity needs to be
added with the ``force_unlazy`` option.)

Two functions are necessary.  One will operate on arrays of data, either fed
from each grid individually or fed from the entire data object at once.  The
second one takes the results of the first, either as lists of arrays or as
single arrays, and returns the final values.  For an example, we look at the
``TotalMass`` function:

.. code-block:: python

   def _TotalMass(data):
       baryon_mass = data["cell_mass"].sum()
       particle_mass = data["ParticleMassMsun"].sum()
       return baryon_mass, particle_mass
   def _combTotalMass(data, baryon_mass, particle_mass):
       return baryon_mass.sum() + particle_mass.sum()
   add_quantity("TotalMass", function=_TotalMass,
                combine_function=_combTotalMass, n_ret = 2)

Once the two functions have been defined, we then call :func:`add_quantity` to
tell it the function that defines the data, the collator function, and the
number of values that get passed between them.  In this case we return both the
particle and the baryon mass, so we have two total values passed from the main
function into the collator.

.. _field_cuts:

Cutting Objects by Field Values
-------------------------------

Data objects can be cut by their field values using the ``cut_region`` 
method.  For example, this could be used to compute the total gas mass within
a certain temperature range, as in the following example.

.. notebook-cell::

   import yt
   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   ad = ds.all_data()
   total_mass = ad.quantities.total_quantity('cell_mass')
   # now select only gas with 1e5 K < T < 1e7 K.
   new_region = ad.cut_region(['obj["temperature"] > 1e5',
                               'obj["temperature"] < 1e7'])
   cut_mass = new_region.quantities.total_quantity('cell_mass')
   print "The fraction of mass in this temperature range is %f." % \
     (cut_mass / total_mass)

The ``cut_region`` function generates a new object containing only the cells 
that meet all of the specified criteria.  The sole argument to ``cut_region`` 
is a list of strings, where each string is evaluated with an ``eval`` 
statement.  ``eval`` is a native Python function that evaluates a string as 
a Python expression.  Any type of data object can be cut with ``cut_region``.  
Objects generated with ``cut_region`` can be used in the same way as all 
other data objects.  For example, a cut region can be visualized by giving 
it as a data_source to a projection.

.. python-script::

   import yt
   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   ad = ds.all_data()
   new_region = ad.cut_region(['obj["density"] > 1e-29'])
   plot = yt.ProjectionPlot(ds, "x", "density", weight_field="density",
                            data_source=new_region)
   plot.save()

.. _extracting-connected-sets:

Connected Sets
--------------

The underlying machinery used in :ref:`clump_finding` is accessible from any
data object.  This includes the ability to obtain and examine topologically
connected sets.  These sets are identified by examining cells between two
threshold values and connecting them.  What is returned to the user is a list
of the intervals of values found, and extracted regions that contain only those
cells that are connected.

To use this, call
:meth:`~yt.data_objects.data_containers.AMR3DData.extract_connected_sets` on
any 3D data object.  This requests a field, the number of levels of levels sets to
extract, the min and the max value between which sets will be identified, and
whether or not to conduct it in log space.

.. code-block:: python

   sp = ds.sphere("max", (1.0, 'pc'))
   contour_values, connected_sets = sp.extract_connected_sets(
        "density", 3, 1e-30, 1e-20)

The first item, ``contour_values``, will be an array of the min value for each
set of level sets.  The second (``connected_sets``) will be a dict of dicts.
The key for the first (outer) dict is the level of the contour, corresponding
to ``contour_values``.  The inner dict returned is keyed by the contour ID.  It
contains :class:`~yt.data_objects.data_containers.AMRExtractedRegionBase`
objects.  These can be queried just as any other data object.

.. _extracting-isocontour-information:

Extracting Isocontour Information
---------------------------------

``yt`` contains an implementation of the `Marching Cubes
<http://en.wikipedia.org/wiki/Marching_cubes>`_ algorithm, which can operate on
3D data objects.  This provides two things.  The first is to identify
isocontours and return either the geometry of those isocontours or to return
another field value sampled along that isocontour.  The second piece of
functionality is to calculate the flux of a field over an isocontour.

Note that these isocontours are not guaranteed to be topologically connected.
In fact, inside a given data object, the marching cubes algorithm will return
all isocontours, not just a single connected one.  This means if you encompass
two clumps of a given density in your data object and extract an isocontour at
that density, it will include both of the clumps.

To extract geometry or sample a field, call
:meth:`~yt.data_objects.data_containers.AMR3DData.extract_isocontours`.  To
calculate a flux, call
:meth:`~yt.data_objects.data_containers.AMR3DData.calculate_isocontour_flux`.
both of these operations will run in parallel.

.. _object-serialization:

Storing and Loading Objects
---------------------------

Often, when operating interactively or via the scripting interface, it is
convenient to save an object or multiple objects out to disk and then restart
the calculation later.  Personally, I found this most useful when dealing with
identification of clumps and contours (see :ref:`cookbook` for a recipe on how
to find clumps and the API documentation for both 
:mod:`~yt.analysis_modules.level_sets.contour_finder.identify_contours`
and :mod:`~yt.analysis_modules.level_sets.clump_handling.Clump`) where 
the identification step can be quite time-consuming, but the analysis 
may be relatively fast.

Typically, the save and load operations are used on 3D data objects.  ``yt``
has a separate set of serialization operations for 2D objects such as
projections.

.. _parameter_file_serialization:

``yt`` will save out 3D objects to disk under the presupposition that the
construction of the objects is the difficult part, rather than the generation
of the data -- this means that you can save out an object as a description of
how to recreate it in space, but not the actual data arrays affiliated with
that object.  The information that is saved includes the dataset off of
which the object "hangs."  It is this piece of information that is the most
difficult; the object, when reloaded, must be able to reconstruct a parameter
file from whatever limited information it has in the save file.

To do this, ``yt`` is able to identify datasets based on a "hash"
generated from the base file name, the "CurrentTimeIdentifier", and the
simulation time.  These three characteristics should never be changed outside
of a simulation, they are independent of the file location on disk, and in
conjunction they should be uniquely identifying.  (This process is all done in
:mod:`~yt.utilities.ParameterFileStorage` via :class:`~yt.utilities.ParameterFileStorage.ParameterFileStore`.)

To save an object, you can either save it in the ``.yt`` file affiliated with
the index or as a standalone file.  For instance, using
:meth:`~yt.data_objects.index.save_object` we can save a sphere.

.. code-block:: python

   import yt
   ds = yt.load("my_data")
   sp = ds.sphere([0.5, 0.5, 0.5], 10.0/ds['kpc'])

   ds.save_object(sp, "sphere_to_analyze_later")


In a later session, we can load it using
:meth:`~yt.data_objects.index.load_object`:

.. code-block:: python

   import yt

   ds = yt.load("my_data")
   sphere_to_analyze = ds.load_object("sphere_to_analyze_later")

Additionally, if we want to store the object independent of the ``.yt`` file,
we can save the object directly:

.. code-block:: python

   import yt

   ds = yt.load("my_data")
   sp = ds.sphere([0.5, 0.5, 0.5], 10.0/ds['kpc'])

   sp.save_object("my_sphere", "my_storage_file.cpkl")

This will store the object as ``my_sphere`` in the file
``my_storage_file.cpkl``, which will be created or accessed using the standard
python module :mod:`shelve`.  Note that if a filename is not supplied, it will
be saved via the index, as above.

To re-load an object saved this way, you can use the shelve module directly:

.. code-block:: python

   import yt
   import shelve

   ds = yt.load("my_data") # not necessary if storeparameterfiles is on

   obj_file = shelve.open("my_storage_file.cpkl")
   ds, obj = obj_file["my_sphere"]

If you have turned on ``storeparameterfiles`` in your configuration,
you won't need to load the parameterfile again, as the load process
will actually do that for you in that case.  Additionally, we can
store multiple objects in a single shelve file, so we have to call the
sphere by name.

.. note:: It's also possible to use the standard :mod:`cPickle` module for
          loading and storing objects -- so in theory you could even save a
          list of objects!

This method works for clumps, as well, and the entire clump index will be
stored and restored upon load.

.. _accessing-fields:

Accessing Fields in Objects
---------------------------

``yt`` utilizes load-on-demand objects to represent physical regions in space.
(see :ref:`how-yt-thinks-about-data`.)  Data objects in ``yt`` all respect the following
protocol for accessing data:

.. code-block:: python

   my_object["density"]

where ``"density"`` can be any field name and ``"my_object"`` any one of
the possible data containers listed at :ref:`available-objects`. For
example, if we wanted to look at the temperature of cells within a
spherical region of radius 10 kpc, centered at [0.5, 0.5, 0.5] in our
simulation box, we would create a sphere object with:

.. code-block:: python

   sp = ds.sphere([0.5, 0.5, 0.5], 10.0/ds['kpc'])

and then look at the temperature of its cells within it via:

.. code-block:: python

   print sp["temperature"]

Information about how to create a new type of object can be found in
:ref:`creating-objects`. The field is returned as a single, flattened
array without spatial information.  The best mechanism for
manipulating spatial data is the :class:`~yt.data_objects.data_containers.AMRCoveringGridBase` object.

The full list of fields that are available can be found as a property of the
Hierarchy or Static Output object that you wish to access.  This property is
calculated every time the object is instantiated.  The full list of fields that
have been identified in the output file, which need no processing (besides unit
conversion) are in the property ``field_list`` and the full list of
potentially-accessible derived fields is available in the property
``derived_field_list``.  You can see these by examining the two properties:

.. code-block:: python

   ds = yt.load("my_data")
   print ds.field_list
   print ds.derived_field_list

When a field is added, it is added to a container that hangs off of the
dataset, as well.  All of the field creation options
(:ref:`derived-field-options`) are accessible through this object:

.. code-block:: python

   ds = yt.load("my_data")
   print ds.field_info["pressure"].get_units()

This is a fast way to examine the units of a given field, and additionally you
can use :meth:`yt.utilities.pydot.get_source` to get the source code:

.. code-block:: python

   field = ds.field_info["pressure"]
   print field.get_source()



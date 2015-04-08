.. _data-objects:

Data Objects
============

What are Data Objects in yt?
----------------------------

Data objects (also called *Data Containers*) are used in yt as convenience 
structures for grouping data in logical ways that make sense in the context 
of the dataset as a whole.  Some of the data objects are geometrical groupings 
of data (e.g. sphere, box, cylinder, etc.).  Others represent 
data products derived from your dataset (e.g. slices, streamlines, surfaces).
Still other data objects group multiple objects together or filter them
(e.g. data collection, cut region).  

To generate standard plots, objects rarely need to be directly constructed.
However, for detailed data inspection as well as hand-crafted derived data,
objects can be exceptionally useful and even necessary.

How to Create and Use an Object
-------------------------------

To create an object, you usually only need a loaded dataset, the name of 
the object type, and the relevant parameters for your object.  Here is a common
example for creating a ``Region`` object that covers all of your data volume.

.. code-block:: python

   import yt
   ds = yt.load("RedshiftOutput0005")
   ad = ds.all_data()

Alternatively, we could create a sphere object of radius 1 kpc on location 
[0.5, 0.5, 0.5]:

.. code-block:: python

   import yt
   ds = yt.load("RedshiftOutput0005")
   sp = ds.sphere([0.5, 0.5, 0.5], (1, 'kpc'))

After an object has been created, it can be used as a data_source to certain
tasks like ``ProjectionPlot`` (see 
:class:`~yt.visualization.plot_window.ProjectionPlot`), one can compute the 
bulk quantities associated with that object (see :ref:`derived-quantities`), 
or the data can be examined directly. For example, if you want to figure out 
the temperature at all indexed locations in the central sphere of your 
dataset you could:

.. code-block:: python

   import yt
   ds = yt.load("RedshiftOutput0005")
   sp = ds.sphere([0.5, 0.5, 0.5], (1, 'kpc'))

   # Show all temperature values
   print sp["temperature"]

   # Print things in a more human-friendly manner: one temperature at a time
   print "(x,  y,  z) Temperature"
   print "-----------------------"
   for i in range(sp["temperature"].size):
       print "(%f,  %f,  %f)    %f" % (sp["x"][i], sp["y"][i], sp["z"][i], sp["temperature"][i])

.. _available-objects:

Available Objects
-----------------

As noted above, there are numerous types of objects.  Here we group them
into:

* *Geometric Objects*
  Data is selected based on spatial shapes in the dataset
* *Filtering Objects*
  Data is selected based on other field criteria
* *Collection Objects*
  Multiple objects grouped together
* *Construction Objects* 
  Objects represent some sort of data product constructed by additional analysis

If you want to create your own custom data object type, see 
:ref:`creating-objects`.

.. _geometric-objects:

Geometric Objects
^^^^^^^^^^^^^^^^^

For 0D, 1D, and 2D geometric objects, if the extent of the object
intersects a grid cell, then the cell is included in the object; however,
for 3D objects the *center* of the cell must be within the object in order
for the grid cell to be incorporated.

0D Objects
""""""""""

**Point** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTPointBase`    
    | Usage: ``point(coord, ds=None, field_parameters=None, data_source=None)``
    | A point defined by a single cell at specified coordinates.

1D Objects
""""""""""

**Ray (Axis-Aligned)** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTOrthoRayBase`
    | Usage: ``ortho_ray(axis, coord, ds=None, field_parameters=None, data_source=None)``
    | A line (of data cells) stretching through the full domain 
      aligned with one of the x,y,z axes.  Defined by an axis and a point
      to be intersected.  Please see this 
      :ref:`note about ray data value ordering <ray-data-ordering>`.

**Ray (Arbitrarily-Aligned)** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTRayBase`
    | Usage: ``ray(start_coord, end_coord, ds=None, field_parameters=None, data_source=None)``
    | A line (of data cells) defined by arbitrary start and end coordinates. 
      Please see this 
      :ref:`note about ray data value ordering <ray-data-ordering>`.

2D Objects
""""""""""

**Slice (Axis-Aligned)** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTSliceBase`
    | Usage: ``slice(axis, coord, center=None, ds=None, field_parameters=None, data_source=None)``
    | A plane normal to one of the axes and intersecting a particular 
      coordinate.

**Slice (Arbitrarily-Aligned)** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTCuttingPlaneBase`
    | Usage: ``cutting(normal, coord, north_vector=None, ds=None, field_parameters=None, data_source=None)``
    | A plane normal to a specified vector and intersecting a particular 
      coordinate.

3D Objects
""""""""""

**All Data** 
    | Function :meth:`~yt.data_objects.static_output.Dataset.all_data`
    | Usage: ``all_data(find_max=False)``
    | ``all_data()`` is a wrapper on the Box Region class which defaults to 
      creating a Region covering the entire dataset domain.  It is effectively 
      ``ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)``.

**Box Region** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTRegionBase`
    | Usage: ``region(center, left_edge, right_edge, fields=None, ds=None, field_parameters=None, data_source=None)``
    | Alternatively: ``box(left_edge, right_edge, fields=None, ds=None, field_parameters=None, data_source=None)``
    | A box-like region aligned with the grid axis orientation.  It is 
      defined by a left_edge, a right_edge, and a center.  The left_edge
      and right_edge are the minimum and maximum bounds in the three axes
      respectively.  The center is arbitrary and must only be contained within
      the left_edge and right_edge.  By using the ``box`` wrapper, the center
      is assumed to be the midpoint between the left and right edges.

**Disk/Cylinder** 
    | Class: :class:`~yt.data_objects.selection_data_containers.YTDiskBase`
    | Usage: ``disk(center, normal, radius, height, fields=None, ds=None, field_parameters=None, data_source=None)``
    | A cylinder defined by a point at the center of one of the circular bases,
      a normal vector to it defining the orientation of the length of the
      cylinder, and radius and height values for the cylinder's dimensions.

**Ellipsoid** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTEllipsoidBase`
    | Usage: ``ellipsoid(center, semi_major_axis_length, semi_medium_axis_length, semi_minor_axis_length, semi_major_vector, tilt, fields=None, ds=None, field_parameters=None, data_source=None)``
    | An ellipsoid with axis magnitudes set by semi_major_axis_length, 
     semi_medium_axis_length, and semi_minor_axis_length.  semi_major_vector 
     sets the direction of the semi_major_axis.  tilt defines the orientation 
     of the semi-medium and semi_minor axes.

**Sphere** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTSphereBase`
    | Usage: ``sphere(center, radius, ds=None, field_parameters=None, data_source=None)``
    | A sphere defined by a central coordinate and a radius.

.. _collection-objects:

Filtering and Collection Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See also the section on :ref:`filtering-data`.

**Intersecting Regions**
    | Most Region objects provide a data_source parameter, which allows you to subselect
    | one region from another (in the coordinate system of the DataSet). Note, this can
    | easily lead to empty data for non-intersecting regions.
    | Usage: ``slice(axis, coord, ds, data_source=sph)``

**Boolean Regions** 
    | **Note: not yet implemented in yt 3.0**
    | Usage: ``boolean()``
    | See :ref:`boolean_data_objects`.

**Filter** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTCutRegionBase`
    | Usage: ``cut_region(base_object, conditionals, ds=None, field_parameters=None)``
    | A ``cut_region`` is a filter which can be applied to any other data 
      object.  The filter is defined by the conditionals present, which 
      apply cuts to the data in the object.  A ``cut_region`` will work
      for either particle fields or mesh fields, but not on both simulaneously.
      For more detailed information and examples, see :ref:`cut-regions`.

**Collection of Data Objects** 
    | Class :class:`~yt.data_objects.selection_data_containers.YTDataCollectionBase`
    | Usage: ``data_collection(center, obj_list, ds=None, field_parameters=None)``
    | A ``data_collection`` is a list of data objects that can be 
      sampled and processed as a whole in a single data object.

.. _construction-objects:

Construction Objects
^^^^^^^^^^^^^^^^^^^^

**Fixed-Resolution Region** 
    | Class :class:`~yt.data_objects.construction_data_containers.YTCoveringGridBase`
    | Usage: ``covering_grid(level, left_edge, dimensions, fields=None, ds=None, num_ghost_zones=0, use_pbar=True, field_parameters=None)``
    | A 3D region with all data extracted to a single, specified resolution.
      See :ref:`examining-grid-data-in-a-fixed-resolution-array`.

**Fixed-Resolution Region with Smoothing** 
    | Class :class:`~yt.data_objects.construction_data_containers.YTSmoothedCoveringGridBase`
    | Usage: ``smoothed_covering_grid(level, left_edge, dimensions, fields=None, ds=None, num_ghost_zones=0, use_pbar=True, field_parameters=None)``
    | A 3D region with all data extracted and interpolated to a single, 
      specified resolution.  Identical to covering_grid, except that it 
      interpolates as necessary from coarse regions to fine.  See 
      :ref:`examining-grid-data-in-a-fixed-resolution-array`.

**Fixed-Resolution Region for Particle Deposition** 
    | Class :class:`~yt.data_objects.construction_data_containers.YTArbitraryGridBase`
    | Usage: ``arbitrary_grid(left_edge, right_edge, dimensions, ds=None, field_parameters=None)``
    | When particles are deposited on to mesh fields, they use the existing
      mesh structure, but this may have too much or too little resolution
      relative to the particle locations (or it may not exist at all!).  An
      `arbitrary_grid` provides a means for generating a new independent mesh 
      structure for particle deposition.  See :ref:`arbitrary-grid` for more 
      information.

**Projection** 
    | Class :class:`~yt.data_objects.construction_data_containers.YTQuadTreeProjBase`
    | Usage: ``proj(field, axis, weight_field=None, center=None, ds=None, data_source=None, method="integrate", field_parameters=None)``
    | A 2D projection of a 3D volume along one of the axis directions.  
      By default, this is a line integral through the entire simulation volume 
      (although it can be a subset of that volume specified by a data object
      with the ``data_source`` keyword).  Alternatively, one can specify 
      a weight_field and different ``method`` values to change the nature
      of the projection outcome.  See :ref:`projection-types` for more information.

**Streamline** 
    | Class :class:`~yt.data_objects.construction_data_containers.YTStreamlineBase`
    | Usage: ``streamline(coord_list, length, fields=None, ds=None, field_parameters=None)``
    | A ``streamline`` can be traced out by identifying a starting coordinate (or 
      list of coordinates) and allowing it to trace a vector field, like gas
      velocity.  See :ref:`streamlines` for more information.

**Surface** 
    | Class :class:`~yt.data_objects.construction_data_containers.YTSurfaceBase`
    | Usage: ``surface(data_source, field, field_value)``
    | The surface defined by all an isocontour in any mesh field.  An existing 
      data object must be provided as the source, as well as a mesh field
      and the value of the field which you desire the isocontour.  See 
      :ref:`extracting-isocontour-information`.

.. _derived-quantities:

Processing Objects: Derived Quantities
--------------------------------------

Derived quantities are a way of calculating some bulk quantities associated
with all of the grid cells contained in a data object.  
Derived quantities can be accessed via the ``quantities`` interface.
Here is an example of how to get the angular momentum vector calculated from 
all the cells contained in a sphere at the center of our dataset.

.. code-block:: python

   ds = load("my_data")
   sp = ds.sphere('c', (10, 'kpc'))
   print sp.quantities.angular_momentum_vector()

Available Derived Quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Angular Momentum Vector**
    | Class :class:`~yt.data_objects.derived_quantities.AngularMomentumVector`
    | Usage: ``angular_momentum_vector(use_gas=True, use_particles=True)``
    | The mass-weighted average angular momentum vector of the particles, gas, 
      or both.

**Bulk Velocity**
    | Class :class:`~yt.data_objects.derived_quantities.BulkVelocity`
    | Usage: ``bulk_velocity(use_gas=True, use_particles=True)``
    | The mass-weighted average velocity of the particles, gas, or both.

**Center of Mass**
    | Class :class:`~yt.data_objects.derived_quantities.CenterOfMass`
    | Usage: ``center_of_mass(use_cells=True, use_particles=False)``
    | The location of the center of mass. By default, it computes of 
      the *non-particle* data in the object, but it can be used on 
      particles, gas, or both.

**Extrema**
    | Class :class:`~yt.data_objects.derived_quantities.Extrema`
    | Usage: ``extrema(fields, non_zero=False)``
    | The extrema of a field or list of fields.

**Maximum Location**
    | Class :class:`~yt.data_objects.derived_quantities.MaxLocation`
    | Usage: ``max_location(fields)``
    | The maximum of a field or list of fields as well
      as the x,y,z location of that maximum.

**Minimum Location**
    | Class :class:`~yt.data_objects.derived_quantities.MinLocation`
    | Usage: ``min_location(fields)``
    | The minimum of a field or list of fields as well
      as the x,y,z location of that minimum.

**Spin Parameter**
    | Class :class:`~yt.data_objects.derived_quantities.SpinParameter`
    | Usage: ``spin_parameter(use_gas=True, use_particles=True)``
    | The spin parameter for the baryons using the particles, gas, or both.

**Total Mass**
    | Class :class:`~yt.data_objects.derived_quantities.TotalMass`
    | Usage: ``total_mass()``
    | The total mass of the object as a tuple of (total gas, total particle)
      mass.

**Total of a Field**
    | Class :class:`~yt.data_objects.derived_quantities.TotalQuantity`
    | Usage: ``total_quantity(fields)``
    | The sum of a given field (or list of fields) over the entire object.

**Weighted Average of a Field**
    | Class :class:`~yt.data_objects.derived_quantities.WeightedAverageQuantity`
    | Usage: ``weighted_average_quantity(fields, weight)``
    | The weighted average of a field (or list of fields)
      over an entire data object.  If you want an unweighted average, 
      then set your weight to be the field: ``ones``.

**Weighted Variance of a Field**
    | Class :class:`~yt.data_objects.derived_quantities.WeightedVariance`
    | Usage: ``weighted_variance(fields, weight)``
    | The weighted variance of a field (or list of fields)
      over an entire data object and the weighted mean.  
      If you want an unweighted variance, then 
      set your weight to be the field: ``ones``.

.. _arbitrary-grid:

Arbitrary Grids Objects for Particle Deposition
-----------------------------------------------

The covering grid and smoothed covering grid objects mandate that they be
exactly aligned with the mesh.  This is a
holdover from the time when yt was used exclusively for data that came in
regularly structured grid patches, and does not necessarily work as well for
data that is composed of discrete objects like particles.  To augment this, the
:class:`~yt.data_objects.construction_data_containers.YTArbitraryGridBase` object 
was created, which enables construction of meshes (onto which particles can be
deposited or smoothed) in arbitrary regions.  This eliminates any assumptions
on yt's part about how the data is organized, and will allow for more
fine-grained control over visualizations.

An example of creating an arbitrary grid would be to construct one, then query
the deposited particle density, like so:

.. code-block:: python

   import yt
   ds = yt.load("snapshot_010.hdf5")

   obj = ds.arbitrary_grid([0.0, 0.0, 0.0], [0.99, 0.99, 0.99],
                          dims=[128, 128, 128])
   print obj["deposit", "all_density"]

While these cannot yet be used as input to projections or slices, slices and
projections can be taken of the data in them and visualized by hand.

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

.. _extracting-connected-sets:

Connected Sets and Clump Finding
--------------------------------

The underlying machinery used in :ref:`clump_finding` is accessible from any
data object.  This includes the ability to obtain and examine topologically
connected sets.  These sets are identified by examining cells between two
threshold values and connecting them.  What is returned to the user is a list
of the intervals of values found, and extracted regions that contain only those
cells that are connected.

To use this, call
:meth:`~yt.data_objects.data_containers.YTSelectionContainer3D.extract_connected_sets` on
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
contains :class:`~yt.data_objects.selection_data_containers.YTCutRegionBase`
objects.  These can be queried just as any other data object.  The clump finder 
(:ref:`clump_finding`) differs from the above method in that the contour 
identification is performed recursively within each individual structure, and 
structures can be kept or remerged later based on additional criteria, such as 
gravitational boundedness.

.. _object-serialization:

Storing and Loading Objects
---------------------------

Often, when operating interactively or via the scripting interface, it is
convenient to save an object or multiple objects out to disk and then restart
the calculation later.  For example, this is useful after clump finding 
(:ref:`clump_finding`), which can be very time consuming.  
Typically, the save and load operations are used on 3D data objects.  yt
has a separate set of serialization operations for 2D objects such as
projections.

yt will save out objects to disk under the presupposition that the
construction of the objects is the difficult part, rather than the generation
of the data -- this means that you can save out an object as a description of
how to recreate it in space, but not the actual data arrays affiliated with
that object.  The information that is saved includes the dataset off of
which the object "hangs."  It is this piece of information that is the most
difficult; the object, when reloaded, must be able to reconstruct a dataset
from whatever limited information it has in the save file.

You can save objects to an output file using the function 
:func:`~yt.data_objects.index.save_object`: 

.. code-block:: python

   import yt
   ds = yt.load("my_data")
   sp = ds.sphere([0.5, 0.5, 0.5], (10.0, 'kpc'))
   sp.save_object("sphere_name", "save_file.cpkl")

This will store the object as ``sphere_name`` in the file
``save_file.cpkl``, which will be created or accessed using the standard
python module :mod:`shelve`.  

To re-load an object saved this way, you can use the shelve module directly:

.. code-block:: python

   import yt
   import shelve
   ds = yt.load("my_data") 
   saved_fn = shelve.open("save_file.cpkl")
   ds, sp = saved_fn["sphere_name"]

Additionally, we can store multiple objects in a single shelve file, so we 
have to call the sphere by name.

For certain data objects such as projections, serialization can be performed
automatically if ``serialize`` option is set to ``True`` in
:ref:`configuration-file`: or set directly in the script:

.. code-block:: python

   from yt.config import ytcfg; ytcfg["yt", "serialize"] = "True"

.. note:: Use serialization with caution. Enabling serialization means that
   once a projection of a dataset has been created (and stored in the .yt file
   in the same directory), any subsequent changes to that dataset will be
   ignored when attempting to create the same projection. So if you take a
   density projection of your dataset in the 'x' direction, then somehow tweak
   that dataset significantly, and take the density projection again, yt will
   default to finding the original projection and 
   :ref:`not your new one <faq-old-data>`:.

.. note:: It's also possible to use the standard :mod:`cPickle` module for
          loading and storing objects -- so in theory you could even save a
          list of objects!

This method works for clumps, as well, and the entire clump index will be
stored and restored upon load.

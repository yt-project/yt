.. _Data-objects:

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
   sp = ds.sphere([0.5, 0.5, 0.5], (1, "kpc"))

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
   sp = ds.sphere([0.5, 0.5, 0.5], (1, "kpc"))

   # Show all temperature values
   print(sp["gas", "temperature"])

   # Print things in a more human-friendly manner: one temperature at a time
   print("(x,  y,  z) Temperature")
   print("-----------------------")
   for i in range(sp["gas", "temperature"].size):
       print(
           "(%f,  %f,  %f)    %f"
           % (
               sp["gas", "x"][i],
               sp["gas", "y"][i],
               sp["gas", "z"][i],
               sp["gas", "temperature"][i],
           )
       )

Data objects can also be cloned; for instance:

.. code-block:: python

   import yt

   ds = yt.load("RedshiftOutput0005")
   sp = ds.sphere([0.5, 0.5, 0.5], (1, "kpc"))
   sp_copy = sp.clone()

This can be useful for when manually chunking data or exploring different field
parameters.

.. _quickly-selecting-data:

Slicing Syntax for Selecting Data
---------------------------------

yt provides a mechanism for easily selecting data while doing interactive work
on the command line.  This allows for region selection based on the full domain
of the object.  Selecting in this manner is exposed through a slice-like
syntax.  All of these attributes are exposed through the ``RegionExpression``
object, which is an attribute of a ``DataSet`` object, called ``r``.

Getting All The Data
^^^^^^^^^^^^^^^^^^^^

The ``.r`` attribute serves as a persistent means of accessing the full data
from a dataset.  You can access this shorthand operation by querying any field
on the ``.r`` object, like so:

.. code-block:: python

   ds = yt.load("RedshiftOutput0005")
   rho = ds.r["gas", "density"]

This will return a *flattened* array of data.  The region expression object
(``r``) doesn't have any derived quantities on it.  This is completely
equivalent to this set of statements:

.. code-block:: python

   ds = yt.load("RedshiftOutput0005")
   dd = ds.all_data()
   rho = dd["gas", "density"]

.. warning::

   One thing to keep in mind with accessing data in this way is that it is
   *persistent*.  It is loaded into memory, and then retained until the dataset
   is deleted or garbage collected.

Selecting Multiresolution Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To select rectilinear regions, where the data is selected the same way that it
is selected in a :ref:`region-reference`, you can utilize slice-like syntax,
supplying start and stop, but not supplying a step argument.  This requires
that three components of the slice must be specified.  These take a start and a
stop, and are for the three axes in simulation order (if your data is ordered
z, y, x for instance, this would be in z, y, x order).

The slices can have both position and, optionally, unit values.  These define
the value with respect to the ``domain_left_edge`` of the dataset.  So for
instance, you could specify it like so:

.. code-block:: python

   ds.r[(100, "kpc"):(200, "kpc"), :, :]

This would return a region that included everything between 100 kpc from the
left edge of the dataset to 200 kpc from the left edge of the dataset in the
first dimension, and which spans the entire dataset in the second and third
dimensions.  By default, if the units are unspecified, they are in the "native"
code units of the dataset.

This works in all types of datasets, as well.  For instance, if you have a
geographic dataset (which is usually ordered latitude, longitude, altitude) you
can easily select, for instance, one hemisphere with a region selection:

.. code-block:: python

   ds.r[:, -180:0, :]

If you specify a single slice, it will be repeated along all three dimensions.
For instance, this will give all data:

.. code-block:: python

   ds.r[:]

And this will select a box running from 0.4 to 0.6 along all three
dimensions:

.. code-block:: python

   ds.r[0.4:0.6]


.. _arbitrary-grid-selection:

Selecting Fixed Resolution Regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yt also provides functionality for selecting regions that have been turned into
voxels.  This returns an :ref:`arbitrary-grid` object.  It can be created by
specifying a complex slice "step", where the start and stop follow the same
rules as above.  This is similar to how the numpy ``mgrid`` operation works.
For instance, this code block will generate a grid covering the full domain,
but converted to being 21x35x100 dimensions:

.. code-block:: python

   region = ds.r[::21j, ::35j, ::100j]

The left and right edges, as above, can be specified to provide bounds as well.
For instance, to select a 10 meter cube, with 24 cells in each dimension, we
could supply:

.. code-block:: python

   region = ds.r[(20, "m"):(30, "m"):24j, (30, "m"):(40, "m"):24j, (7, "m"):(17, "m"):24j]

This can select both particles and mesh fields.  Mesh fields will be 3D arrays,
and generated through volume-weighted overlap calculations.

Selecting Slices
^^^^^^^^^^^^^^^^

If one dimension is specified as a single value, that will be the dimension
along which a slice is made.  This provides a simple means of generating a
slice from a subset of the data.  For instance, to create a slice of a dataset,
you can very simply specify the full domain along two axes:

.. code-block:: python

    sl = ds.r[:, :, 0.25]

This can also be very easily plotted:

.. code-block:: python

   sl = ds.r[:, :, 0.25]
   sl.plot()

This accepts arguments the same way:

.. code-block:: python

   sl = ds.r[(20.1, "km"):(31.0, "km"), (504.143, "m"):(1000.0, "m"), (900.1, "m")]
   sl.plot()

Making Image Buffers
^^^^^^^^^^^^^^^^^^^^

Using the slicing syntax above for choosing a slice, if you also provide an
imaginary step value you can obtain a
:class:`~yt.visualization.api.FixedResolutionBuffer` of the chosen resolution.

For instance, to obtain a 1024 by 1024 buffer covering the entire
domain but centered at 0.5 in code units, you can do:

.. code-block:: python

   frb = ds.r[0.5, ::1024j, ::1024j]

This ``frb`` object then can be queried like a normal fixed resolution buffer,
and it will return arrays of shape (1024, 1024).

Making Rays
^^^^^^^^^^^

The slicing syntax can also be used select 1D rays of points, whether along
an axis or off-axis. To create a ray along an axis:

.. code-block:: python

   ortho_ray = ds.r[(500.0, "kpc"), (200, "kpc"):(300.0, "kpc"), (-2.0, "Mpc")]

To create a ray off-axis, use a single slice between the start and end points
of the ray:

.. code-block:: python

   start = [0.1, 0.2, 0.3]  # interpreted in code_length
   end = [0.4, 0.5, 0.6]  # interpreted in code_length
   ray = ds.r[start:end]

As for the other slicing options, combinations of unitful quantities with even
different units can be used. Here's a somewhat convoluted (yet working) example:

.. code-block:: python

   start = ((500.0, "kpc"), (0.2, "Mpc"), (100.0, "kpc"))
   end = ((1.0, "Mpc"), (300.0, "kpc"), (0.0, "kpc"))
   ray = ds.r[start:end]

Making Fixed-Resolution Rays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rays can also be constructed to have fixed resolution if an imaginary step value
is provided, similar to the 2 and 3-dimensional cases described above. This
works for rays directed along an axis:

.. code-block:: python

   ortho_ray = ds.r[0.1:0.6:500j, 0.3, 0.2]

or off-axis rays as well:

.. code-block:: python

   start = [0.1, 0.2, 0.3]  # interpreted in code_length
   end = [0.4, 0.5, 0.6]  # interpreted in code_length
   ray = ds.r[start:end:100j]

Selecting Points
^^^^^^^^^^^^^^^^

Finally, you can quickly select a single point within the domain by providing
a single coordinate for every axis:

.. code-block:: python

   pt = ds.r[(10.0, "km"), (200, "m"), (1.0, "km")]

Querying this object for fields will give you the value of the field at that
point.

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
    | Class :class:`~yt.data_objects.selection_data_containers.YTPoint`
    | Usage: ``point(coord, ds=None, field_parameters=None, data_source=None)``
    | A point defined by a single cell at specified coordinates.

1D Objects
""""""""""

**Ray (Axis-Aligned)**
    | Class :class:`~yt.data_objects.selection_data_containers.YTOrthoRay`
    | Usage: ``ortho_ray(axis, coord, ds=None, field_parameters=None, data_source=None)``
    | A line (of data cells) stretching through the full domain
      aligned with one of the x,y,z axes.  Defined by an axis and a point
      to be intersected.  Please see this
      :ref:`note about ray data value ordering <ray-data-ordering>`.

**Ray (Arbitrarily-Aligned)**
    | Class :class:`~yt.data_objects.selection_data_containers.YTRay`
    | Usage: ``ray(start_coord, end_coord, ds=None, field_parameters=None, data_source=None)``
    | A line (of data cells) defined by arbitrary start and end coordinates.
      Please see this
      :ref:`note about ray data value ordering <ray-data-ordering>`.

2D Objects
""""""""""

**Slice (Axis-Aligned)**
    | Class :class:`~yt.data_objects.selection_data_containers.YTSlice`
    | Usage: ``slice(axis, coord, center=None, ds=None, field_parameters=None, data_source=None)``
    | A plane normal to one of the axes and intersecting a particular
      coordinate.

**Slice (Arbitrarily-Aligned)**
    | Class :class:`~yt.data_objects.selection_data_containers.YTCuttingPlane`
    | Usage: ``cutting(normal, coord, north_vector=None, ds=None, field_parameters=None, data_source=None)``
    | A plane normal to a specified vector and intersecting a particular
      coordinate.

.. _region-reference:

3D Objects
""""""""""

**All Data**
    | Function :meth:`~yt.data_objects.static_output.Dataset.all_data`
    | Usage: ``all_data(find_max=False)``
    | ``all_data()`` is a wrapper on the Box Region class which defaults to
      creating a Region covering the entire dataset domain.  It is effectively
      ``ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)``.

**Box Region**
    | Class :class:`~yt.data_objects.selection_data_containers.YTRegion`
    | Usage: ``region(center, left_edge, right_edge, fields=None, ds=None, field_parameters=None, data_source=None)``
    | Alternatively: ``box(left_edge, right_edge, fields=None, ds=None, field_parameters=None, data_source=None)``
    | A box-like region aligned with the grid axis orientation.  It is
      defined by a left_edge, a right_edge, and a center.  The left_edge
      and right_edge are the minimum and maximum bounds in the three axes
      respectively.  The center is arbitrary and must only be contained within
      the left_edge and right_edge.  By using the ``box`` wrapper, the center
      is assumed to be the midpoint between the left and right edges.

**Disk/Cylinder**
    | Class: :class:`~yt.data_objects.selection_data_containers.YTDisk`
    | Usage: ``disk(center, normal, radius, height, fields=None, ds=None, field_parameters=None, data_source=None)``
    | A cylinder defined by a point at the center of one of the circular bases,
      a normal vector to it defining the orientation of the length of the
      cylinder, and radius and height values for the cylinder's dimensions.
      Note: ``height`` is the distance from midplane to the top or bottom of the
      cylinder, i.e., ``height`` is half that of the cylinder object that is
      created.

**Ellipsoid**
    | Class :class:`~yt.data_objects.selection_data_containers.YTEllipsoid`
    | Usage: ``ellipsoid(center, semi_major_axis_length, semi_medium_axis_length, semi_minor_axis_length, semi_major_vector, tilt, fields=None, ds=None, field_parameters=None, data_source=None)``
    | An ellipsoid with axis magnitudes set by ``semi_major_axis_length``,
     ``semi_medium_axis_length``, and ``semi_minor_axis_length``.  ``semi_major_vector``
     sets the direction of the ``semi_major_axis``.  ``tilt`` defines the orientation
     of the semi-medium and semi_minor axes.

**Sphere**
    | Class :class:`~yt.data_objects.selection_data_containers.YTSphere`
    | Usage: ``sphere(center, radius, ds=None, field_parameters=None, data_source=None)``
    | A sphere defined by a central coordinate and a radius.

**Minimal Bounding Sphere**
    | Class :class:`~yt.data_objects.selection_data_containers.YTMinimalSphere`
    | Usage: ``minimal_sphere(points, ds=None, field_parameters=None, data_source=None)``
    | A sphere that contains all the points passed as argument.

.. _collection-objects:

Filtering and Collection Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See also the section on :ref:`filtering-data`.

**Intersecting Regions**
    | Most Region objects provide a data_source parameter, which allows you to subselect
    | one region from another (in the coordinate system of the DataSet). Note, this can
    | easily lead to empty data for non-intersecting regions.
    | Usage: ``slice(axis, coord, ds, data_source=sph)``

**Union Regions**
    | Usage: ``union()``
    | See :ref:`boolean_data_objects`.

**Intersection Regions**
    | Usage: ``intersection()``
    | See :ref:`boolean_data_objects`.

**Filter**
    | Class :class:`~yt.data_objects.selection_data_containers.YTCutRegion`
    | Usage: ``cut_region(base_object, conditionals, ds=None, field_parameters=None)``
    | A ``cut_region`` is a filter which can be applied to any other data
      object.  The filter is defined by the conditionals present, which
      apply cuts to the data in the object.  A ``cut_region`` will work
      for either particle fields or mesh fields, but not on both simultaneously.
      For more detailed information and examples, see :ref:`cut-regions`.

**Collection of Data Objects**
    | Class :class:`~yt.data_objects.selection_data_containers.YTDataCollection`
    | Usage: ``data_collection(center, obj_list, ds=None, field_parameters=None)``
    | A ``data_collection`` is a list of data objects that can be
      sampled and processed as a whole in a single data object.

.. _construction-objects:

Construction Objects
^^^^^^^^^^^^^^^^^^^^

**Fixed-Resolution Region**
    | Class :class:`~yt.data_objects.construction_data_containers.YTCoveringGrid`
    | Usage: ``covering_grid(level, left_edge, dimensions, fields=None, ds=None, num_ghost_zones=0, use_pbar=True, field_parameters=None)``
    | A 3D region with all data extracted to a single, specified resolution.
      See :ref:`examining-grid-data-in-a-fixed-resolution-array`.

**Fixed-Resolution Region with Smoothing**
    | Class :class:`~yt.data_objects.construction_data_containers.YTSmoothedCoveringGrid`
    | Usage: ``smoothed_covering_grid(level, left_edge, dimensions, fields=None, ds=None, num_ghost_zones=0, use_pbar=True, field_parameters=None)``
    | A 3D region with all data extracted and interpolated to a single,
      specified resolution.  Identical to covering_grid, except that it
      interpolates as necessary from coarse regions to fine.  See
      :ref:`examining-grid-data-in-a-fixed-resolution-array`.

**Fixed-Resolution Region**
    | Class :class:`~yt.data_objects.construction_data_containers.YTArbitraryGrid`
    | Usage: ``arbitrary_grid(left_edge, right_edge, dimensions, ds=None, field_parameters=None)``
    | When particles are deposited on to mesh fields, they use the existing
      mesh structure, but this may have too much or too little resolution
      relative to the particle locations (or it may not exist at all!).  An
      `arbitrary_grid` provides a means for generating a new independent mesh
      structure for particle deposition and simple mesh field interpolation.
      See :ref:`arbitrary-grid` for more information.

**Projection**
    | Class :class:`~yt.data_objects.construction_data_containers.YTQuadTreeProj`
    | Usage: ``proj(field, axis, weight_field=None, center=None, ds=None, data_source=None, method="integrate", field_parameters=None)``
    | A 2D projection of a 3D volume along one of the axis directions.
      By default, this is a line integral through the entire simulation volume
      (although it can be a subset of that volume specified by a data object
      with the ``data_source`` keyword).  Alternatively, one can specify
      a weight_field and different ``method`` values to change the nature
      of the projection outcome.  See :ref:`projection-types` for more information.

**Streamline**
    | Class :class:`~yt.data_objects.construction_data_containers.YTStreamline`
    | Usage: ``streamline(coord_list, length, fields=None, ds=None, field_parameters=None)``
    | A ``streamline`` can be traced out by identifying a starting coordinate (or
      list of coordinates) and allowing it to trace a vector field, like gas
      velocity.  See :ref:`streamlines` for more information.

**Surface**
    | Class :class:`~yt.data_objects.construction_data_containers.YTSurface`
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

   import yt

   ds = yt.load("my_data")
   sp = ds.sphere("c", (10, "kpc"))
   print(sp.quantities.angular_momentum_vector())

Some quantities can be calculated for a specific particle type only. For example, to
get the center of mass of only the stars within the sphere:

.. code-block:: python

   import yt

   ds = yt.load("my_data")
   sp = ds.sphere("c", (10, "kpc"))
   print(
       sp.quantities.center_of_mass(
           use_gas=False, use_particles=True, particle_type="star"
       )
   )


Quickly Processing Data
^^^^^^^^^^^^^^^^^^^^^^^

Most data objects now have multiple numpy-like methods that allow you to
quickly process data.  More of these methods will be added over time and added
to this list.  Most, if not all, of these map to other yt operations and are
designed as syntactic sugar to slightly simplify otherwise somewhat obtuse
pipelines.

These operations are parallelized.

You can compute the extrema of a field by using the ``max`` or ``min``
functions.  This will cache the extrema in between, so calling ``min`` right
after ``max`` will be considerably faster.  Here is an example.

.. code-block:: python

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   reg = ds.r[0.3:0.6, 0.2:0.4, 0.9:0.95]
   min_rho = reg.min(("gas", "density"))
   max_rho = reg.max(("gas", "density"))

This is equivalent to:

.. code-block:: python

   min_rho, max_rho = reg.quantities.extrema(("gas", "density"))

The ``max`` operation can also compute the maximum intensity projection:

.. code-block:: python

   proj = reg.max(("gas", "density"), axis="x")
   proj.plot()

This is equivalent to:

.. code-block:: python

   proj = ds.proj(("gas", "density"), "x", data_source=reg, method="mip")
   proj.plot()

The ``min`` operator does not do this, however, as a minimum intensity
projection is not currently implemented.

You can also compute the ``mean`` value, which accepts a field, axis and weight
function.  If the axis is not specified, it will return the average value of
the specified field, weighted by the weight argument.  The weight argument
defaults to ``ones``, which performs an arithmetic average.  For instance:

.. code-block:: python

   mean_rho = reg.mean(("gas", "density"))
   rho_by_vol = reg.mean(("gas", "density"), weight=("gas", "cell_volume"))

This is equivalent to:

.. code-block:: python

   mean_rho = reg.quantities.weighted_average(
       ("gas", "density"), weight_field=("index", "ones")
   )
   rho_by_vol = reg.quantities.weighted_average(
       ("gas", "density"), weight_field=("gas", "cell_volume")
   )

If an axis is provided, it will project along that axis and return it to you:

.. code-block:: python

   rho_proj = reg.mean(("gas", "temperature"), axis="y", weight=("gas", "density"))
   rho_proj.plot()

The ``sum`` function will add all the values in the data object.  It accepts a
field and, optionally, an axis.  If the axis is left unspecified, it will sum
the values in the object:

.. code-block:: python

   vol = reg.sum(("gas", "cell_volume"))

If the axis is specified, it will compute a projection using the method ``sum``
(which does *not* take into account varying path length!) and return that to
you.

.. code-block:: python

   cell_count = reg.sum(("index", "ones"), axis="z")
   cell_count.plot()

To compute a projection where the path length *is* taken into account, you can
use the ``integrate`` function:

.. code-block:: python

   proj = reg.integrate(("gas", "density"), "x")

All of these projections supply the data object as their base input.

Often, it can be useful to sample a field at the minimum and maximum of a
different field.  You can use the ``argmax`` and ``argmin`` operations to do
this.

.. code-block:: python

   reg.argmin(("gas", "density"), axis=("gas", "temperature"))

This will return the temperature at the minimum density.

If you don't specify an ``axis``, it will return the spatial position of
the maximum value of the queried field.  Here is an example::

  x, y, z = reg.argmin(("gas", "density"))

Available Derived Quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Angular Momentum Vector**
    | Class :class:`~yt.data_objects.derived_quantities.AngularMomentumVector`
    | Usage: ``angular_momentum_vector(use_gas=True, use_particles=True, particle_type='all')``
    | The mass-weighted average angular momentum vector of the particles, gas,
      or both. The quantity can be calculated for all particles or a given
      particle_type only.

**Bulk Velocity**
    | Class :class:`~yt.data_objects.derived_quantities.BulkVelocity`
    | Usage: ``bulk_velocity(use_gas=True, use_particles=True, particle_type='all')``
    | The mass-weighted average velocity of the particles, gas, or both.
      The quantity can be calculated for all particles or a given
      particle_type only.

**Center of Mass**
    | Class :class:`~yt.data_objects.derived_quantities.CenterOfMass`
    | Usage: ``center_of_mass(use_cells=True, use_particles=False, particle_type='all')``
    | The location of the center of mass. By default, it computes of
      the *non-particle* data in the object, but it can be used on
      particles, gas, or both. The quantity can be
      calculated for all particles or a given particle_type only.


**Extrema**
    | Class :class:`~yt.data_objects.derived_quantities.Extrema`
    | Usage: ``extrema(fields, non_zero=False)``
    | The extrema of a field or list of fields.

**Maximum Location Sampling**
    | Class :class:`~yt.data_objects.derived_quantities.SampleAtMaxFieldValues`
    | Usage: ``sample_at_max_field_values(fields, sample_fields)``
    | The value of sample_fields at the maximum value in fields.

**Minimum Location Sampling**
    | Class :class:`~yt.data_objects.derived_quantities.SampleAtMinFieldValues`
    | Usage: ``sample_at_min_field_values(fields, sample_fields)``
    | The value of sample_fields at the minimum value in fields.

**Minimum Location**
    | Class :class:`~yt.data_objects.derived_quantities.MinLocation`
    | Usage: ``min_location(fields)``
    | The minimum of a field or list of fields as well
      as the x,y,z location of that minimum.

**Maximum Location**
    | Class :class:`~yt.data_objects.derived_quantities.MaxLocation`
    | Usage: ``max_location(fields)``
    | The maximum of a field or list of fields as well
      as the x,y,z location of that maximum.

**Spin Parameter**
    | Class :class:`~yt.data_objects.derived_quantities.SpinParameter`
    | Usage: ``spin_parameter(use_gas=True, use_particles=True, particle_type='all')``
    | The spin parameter for the baryons using the particles, gas, or both. The
      quantity can be calculated for all particles or a given particle_type only.

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

**Weighted Standard Deviation of a Field**
    | Class :class:`~yt.data_objects.derived_quantities.WeightedStandardDeviation`
    | Usage: ``weighted_standard_deviation(fields, weight)``
    | The weighted standard deviation of a field (or list of fields)
      over an entire data object and the weighted mean.
      If you want an unweighted standard deviation, then
      set your weight to be the field: ``ones``.

.. _arbitrary-grid:

Arbitrary Grids Objects
-----------------------

The covering grid and smoothed covering grid objects mandate that they be
exactly aligned with the mesh.  This is a
holdover from the time when yt was used exclusively for data that came in
regularly structured grid patches, and does not necessarily work as well for
data that is composed of discrete objects like particles.  To augment this, the
:class:`~yt.data_objects.construction_data_containers.YTArbitraryGrid` object
was created, which enables construction of meshes (onto which particles can be
deposited or smoothed) in arbitrary regions.  This eliminates any assumptions
on yt's part about how the data is organized, and will allow for more
fine-grained control over visualizations.

An example of creating an arbitrary grid would be to construct one, then query
the deposited particle density, like so:

.. code-block:: python

   import yt

   ds = yt.load("snapshot_010.hdf5")

   obj = ds.arbitrary_grid([0.0, 0.0, 0.0], [0.99, 0.99, 0.99], dims=[128, 128, 128])
   print(obj["deposit", "all_density"])

While these cannot yet be used as input to projections or slices, slices and
projections can be taken of the data in them and visualized by hand.

These objects, as of yt 3.3, are now also able to "voxelize" mesh fields.  This
means that you can query the "density" field and it will return the density
field as deposited, identically to how it would be deposited in a fixed
resolution buffer.  Note that this means that contributions from misaligned or
partially-overlapping cells are added in a volume-weighted way, which makes it
inappropriate for some types of analysis.

.. _boolean_data_objects:

Combining Objects: Boolean Data Objects
---------------------------------------

A special type of data object is the *boolean* data object, which works with
data selection objects of any dimension.  It is built by relating already existing
data objects with the bitwise operators for AND, OR and XOR, as well as the
subtraction operator.  These are created by using the operators ``&`` for an
intersection ("AND"), ``|`` for a union ("OR"), ``^`` for an exclusive or
("XOR"), and ``+`` and ``-`` for addition ("OR") and subtraction ("NEG").
Here are some examples:

.. code-block:: python

   import yt

   ds = yt.load("snapshot_010.hdf5")

   sp1 = ds.sphere("c", (0.1, "unitary"))
   sp2 = ds.sphere(sp1.center + 2.0 * sp1.radius, (0.2, "unitary"))
   sp3 = ds.sphere("c", (0.05, "unitary"))

   new_obj = sp1 + sp2
   cutout = sp1 - sp3
   sp4 = sp1 ^ sp2
   sp5 = sp1 & sp2


Note that the ``+`` operation and the ``|`` operation are identical.  For when
multiple objects are to be combined in an intersection or a union, there are
the data objects ``intersection`` and ``union`` which can be called, and which
will yield slightly higher performance than a sequence of calls to ``+`` or
``&``.  For instance:

.. code-block:: python

   import yt

   ds = yt.load("Enzo_64/DD0043/data0043")
   sp1 = ds.sphere((0.1, 0.2, 0.3), (0.05, "unitary"))
   sp2 = ds.sphere((0.2, 0.2, 0.3), (0.10, "unitary"))
   sp3 = ds.sphere((0.3, 0.2, 0.3), (0.15, "unitary"))

   isp = ds.intersection([sp1, sp2, sp3])
   usp = ds.union([sp1, sp2, sp3])

The ``isp`` and ``usp`` objects will act the same as a set of chained ``&`` and
``|`` operations (respectively) but are somewhat easier to construct.

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

   sp = ds.sphere("max", (1.0, "pc"))
   contour_values, connected_sets = sp.extract_connected_sets(
       ("gas", "density"), 3, 1e-30, 1e-20
   )

The first item, ``contour_values``, will be an array of the min value for each
set of level sets.  The second (``connected_sets``) will be a dict of dicts.
The key for the first (outer) dict is the level of the contour, corresponding
to ``contour_values``.  The inner dict returned is keyed by the contour ID.  It
contains :class:`~yt.data_objects.selection_data_containers.YTCutRegion`
objects.  These can be queried just as any other data object.  The clump finder
(:ref:`clump_finding`) differs from the above method in that the contour
identification is performed recursively within each individual structure, and
structures can be kept or remerged later based on additional criteria, such as
gravitational boundedness.

.. _object-serialization:

Storing and Loading Objects
---------------------------

Often, when operating interactively or via the scripting interface, it is
convenient to save an object to disk and then restart the calculation later or
transfer the data from a container to another filesystem.  This can be
particularly useful when working with extremely large datasets.  Field data
can be saved to disk in a format that allows for it to be reloaded just like
a regular dataset.  For information on how to do this, see
:ref:`saving-data-containers`.

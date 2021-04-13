.. _low-level-data-inspection:

Low-Level Data Inspection: Accessing Raw Data
=============================================

yt can not only provide high-level access to data, such as through slices,
projections, object queries and the like, but it can also provide low-level
access to the raw data.

.. note:: This section is tuned for patch- or block-based simulations.  Future
          versions of yt will enable more direct access to particle and oct
          based simulations.  For now, these are represented as patches, with
          the attendant properties.

For a more basic introduction, see :ref:`quickstart` and more specifically
:ref:`data_inspection`.

.. _examining-grid-hierarchies:

Examining Grid Hierarchies
--------------------------

yt organizes grids in a hierarchical fashion; a coarser grid that contains (or
overlaps with) a finer grid is referred to as its parent.  yt organizes these
only a single level of refinement at a time.  To access grids, the ``grids``
attribute on a :class:`~yt.geometry.grid_geometry_handler.GridIndex` object.  (For
fast operations, a number of additional arrays prefixed with ``grid`` are also
available, such as ``grid_left_edges`` and so on.)  This returns an instance of
:class:`~yt.data_objects.grid_patch.AMRGridPatch`, which can be queried for
either data or index information.

The :class:`~yt.data_objects.grid_patch.AMRGridPatch` object itself provides
the following attributes:

* ``Children``: a list of grids contained within this one, of one higher level
  of refinement
* ``Parent``: a single object or a list of objects this grid is contained
  within, one level of refinement coarser
* ``child_mask``: a mask of 0's and 1's, representing where no finer data is
  available in refined grids (1) or where this grid is covered by finer regions
  (0).  Note that to get back the final data contained within a grid, one can
  multiple a field by this attribute.
* ``child_indices``: a mask of booleans, where False indicates no finer data
  is available.  This is essentially the inverse of ``child_mask``.
* ``child_index_mask``: a mask of indices into the ``ds.index.grids`` array of the
  child grids.
* ``LeftEdge``: the left edge, in native code coordinates, of this grid
* ``RightEdge``: the right edge, in native code coordinates, of this grid
* ``dds``: the width of a cell in this grid
* ``id``: the id (not necessarily the index) of this grid.  Defined such that
  subtracting the property ``_id_offset`` gives the index into ``ds.index.grids``.
* ``NumberOfParticles``: the number of particles in this grid
* ``OverlappingSiblings``: a list of sibling grids that this grid overlaps
  with.  Likely only defined for Octree-based codes.

In addition, the method
:meth:`~yt.data_objects.grid_patch.AMRGridPatch.get_global_startindex` can be
used to get the integer coordinates of the upper left edge.  These integer
coordinates are defined with respect to the current level; this means that they
are the offset of the left edge, with respect to the left edge of the domain,
divided by the local ``dds``.

To traverse a series of grids, this type of construction can be used:

.. code-block:: python

   g = ds.index.grids[1043]
   g2 = g.Children[1].Children[0]
   print(g2.LeftEdge)

.. _examining-grid-data:

Examining Grid Data
-------------------

Once you have identified a grid you wish to inspect, there are two ways to
examine data.  You can either ask the grid to read the data and pass it to you
as normal, or you can manually intercept the data from the IO handler and
examine it before it has been unit converted.  This allows for much more raw
data inspection.

To access data that has been read in the typical fashion and unit-converted as
normal, you can access the grid as you would a normal object:

.. code-block:: python

   g = ds.index.grids[1043]
   print(g["density"])
   print(g["density"].min())

To access the raw data (as found in the file), use

.. code-block:: python

   g = ds.index.grids[1043]
   rho = g["density"].in_base("code")

.. _finding-data-at-fixed-points:

Finding Data at Fixed Points
----------------------------

One of the most common questions asked of data is, what is the value *at this
specific point*.  While there are several ways to find out the answer to this
question, a few helper routines are provided as well.  To identify the
finest-resolution (i.e., most canonical) data at a given point, use
the point data object::

  from yt.units import kpc
  point_obj = ds.point([30, 75, 80]*kpc)
  density_at_point = point_obj['gas', 'density']

The point data object works just like any other yt data object. It is special
because it is the only zero-dimensional data object: it will only return data at
the exact point specified when creating the point data object. For more
information about yt data objects, see :ref:`Data-objects`.

If you need to find field values at many points, the
:meth:`~yt.data_objects.static_output.Dataset.find_field_values_at_points`
function may be more efficient. This function returns a nested list of field
values at multiple points in the simulation volume. For example, if one wanted
to find the value of a mesh field at the location of the particles in a
simulation, one could do::

  ad = ds.all_data()
  ppos = ad['all', 'particle_position']
  ppos_den_vel = ds.find_field_values_at_points(
      ['density', 'velocity_x'], ppos)

In this example, ``ppos_den_vel`` will be a list of arrays. The first array will
contain the density values at the particle positions, the second will contain
the x velocity values at the particle positions.

.. _examining-grid-data-in-a-fixed-resolution-array:

Examining Grid Data in a Fixed Resolution Array
-----------------------------------------------

If you have a dataset, either AMR or single resolution, and you want to just
stick it into a fixed resolution numpy array for later examination, then you
want to use a :ref:`Covering Grid <available-objects>`.  You must specify the
maximum level at which to sample the data, a left edge of the data where you
will start, and the resolution at which you want to sample.

For example, let's use the :ref:`sample dataset <getting-sample-data>`
``Enzo_64``.  This dataset is at a resolution of 64^3 with 5 levels of AMR,
so if we want a 64^3 array covering the entire volume and sampling just the
lowest level data, we run:

.. code-block:: python

   import yt

   ds = yt.load("Enzo_64/DD0043/data0043")
   all_data_level_0 = ds.covering_grid(level=0, left_edge=[0, 0.0, 0.0], dims=[64, 64, 64])

Note that we can also get the same result and rely on the dataset to know
its own underlying dimensions:

.. code-block:: python

   all_data_level_0 = ds.covering_grid(
       level=0, left_edge=[0, 0.0, 0.0], dims=ds.domain_dimensions
   )

We can now access our underlying data at the lowest level by specifying what
:ref:`field <field-list>` we want to examine:

.. code-block:: python

  print(all_data_level_0["density"].shape)
  # (64, 64, 64)

  print(all_data_level_0["density"])
  # array([[[  1.92588925e-31,   1.74647692e-31,   2.54787518e-31, ...,

  print(all_data_level_0["temperature"].shape)
  # (64, 64, 64)

If you create a covering grid that spans two child grids of a single parent
grid, it will fill those zones covered by a zone of a child grid with the
data from that child grid. Where it is covered only by the parent grid, the
cells from the parent grid will be duplicated (appropriately) to fill the
covering grid.

Let's say we now want to look at that entire data volume and sample it at
a higher resolution (i.e. level 2).  As stated above, we'll be oversampling
under-refined regions, but that's OK.  We must also increase the resolution
of our output array by a factor of 2^2 in each direction to hold this new
larger dataset:

.. code-block:: python

   all_data_level_2 = ds.covering_grid(
       level=2, left_edge=[0, 0.0, 0.0], dims=ds.domain_dimensions * 2 ** 2
   )

And let's see what's the density in the central location:

.. code-block:: python

   print(all_data_level_2["density"].shape)
   (256, 256, 256)

   print(all_data_level_2["density"][128, 128, 128])
   1.7747457571203124e-31

There are two different types of covering grids: unsmoothed and smoothed.
Smoothed grids will be filled through a cascading interpolation process;
they will be filled at level 0, interpolated to level 1, filled at level 1,
interpolated to level 2, filled at level 2, etc. This will help to reduce
edge effects. Unsmoothed covering grids will not be interpolated, but rather
values will be duplicated multiple times.

To sample our dataset from above with a smoothed covering grid in order
to reduce edge effects, it is a nearly identical process:

.. code-block:: python

   all_data_level_2_s = ds.smoothed_covering_grid(
       2, [0.0, 0.0, 0.0], ds.domain_dimensions * 2 ** 2
   )

   print(all_data_level_2_s["density"].shape)
   (256, 256, 256)

   print(all_data_level_2_s["density"][128, 128, 128])
   1.763744852165591e-31

.. _examining-image-data-in-a-fixed-resolution-array:

Examining Image Data in a Fixed Resolution Array
------------------------------------------------

In the same way that one can sample a multi-resolution 3D dataset by placing
it into a fixed resolution 3D array as a
:ref:`Covering Grid <examining-grid-data-in-a-fixed-resolution-array>`, one can
also access the raw image data that is returned from various yt functions
directly as a fixed resolution array.  This provides a means for bypassing the
yt method for generating plots, and allows the user the freedom to use
whatever interface they wish for displaying and saving their image data.
You can use the :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`
to accomplish this as described in :ref:`fixed-resolution-buffers`.

High-level Information about Particles
--------------------------------------

There are a number of high-level helpers attached to ``Dataset`` objects to find
out information about the particles in an output file. First, one can check if
there are any particles in a dataset at all by examining
``ds.particles_exist``. This will be ``True`` for datasets the include particles
and ``False`` otherwise.

One can also see which particle types are available in a dataset. Particle types
that are available in the dataset's on-disk output are known as "raw" particle
types, and they will appear in ``ds.particle_types_raw``. Particle types that
are dynamically defined via a particle filter of a particle union will also
appear in the ``ds.particle_types`` list. If the simulation only has one
particle type on-disk, its name will by ``'io'``. If there is more than one
particle type, the names of the particle types will be inferred from the output
file. For example, Gadget HDF5 files have particle type names like ``PartType0``
and ``PartType1``, while Enzo data, which usually only has one particle type,
will only have a particle named ``io``.

Finally, one can see the number of each particle type by inspecting
``ds.particle_type_counts``. This will be a dictionary mapping the names of
particle types in ``ds.particle_types_raw`` to the number of each particle type
in a simulation output.

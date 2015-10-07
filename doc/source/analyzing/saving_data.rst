.. _saving_data

Saving Reloadable Data
======================

Most of the data loaded into or generated with yt can be saved to a
format that can be reloaded as a first-class dataset.  This includes
the following:

  * geometric data containers (regions, spheres, disks, rays, etc.)

  * grid data containers (covering grids, arbitrary grids, fixed
    resolution buffers)

  * spatial plots (projections, slices, cutting planes)

  * profiles

  * generic array data

In the case of projections, slices, and profiles, reloaded data can be
used to remake plots using the methods decribed in :ref:`how-to-make-plots`.

Geometric Data Containers
-------------------------

Data from geometric data containers can be saved with the
:func:`~yt.data_objects.data_containers.save_as_dataset`` function.

.. code-block:: python

   import yt
   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")

   sphere = ds.sphere([0.5]*3, (10, "Mpc"))
   fn = sphere.save_as_dataset(fields=["density", "particle_mass"])

This function will return the name of the file to which the dataset
was saved.  The filename will be a combination of the name of the
original dataset and the type of data container.  Optionally, a
specific filename can be given with the ``filename`` keyword.  If no
fields are given, the fields that have previously been queried will
be saved.

The newly created dataset can be loaded like all other supported
data through ``yt.load``.  Once loaded, field data can be accessed
through the traditional data containers.  Grid data is accessed by
the ``grid`` data type and particle data is accessed with the
original particle type.  As with the original dataset, grid
positions and cell sizes are accessible with, for example,
("grid", "x") and ("grid", "dx").  Particle positions are
accessible as (<particle_type>, "particle_position_x").  All original
simulation parameters are accessible in the ``parameters``
dictionary, normally associated with all datasets.

.. code-block:: python

   sphere_ds = yt.load("DD0046_sphere.h5")
   ad = sphere_ds.all_data()

   # grid data
   print ad["grid", "density"]
   print ad["grid", "x"]
   print ad["grid", "dx"]

   # particle data
   print ad["all", "particle_mass"]
   print ad["all", "particle_position_x"]

Note that because field data queried from geometric containers is
returned as unordered 1D arrays, data container datasets are treated,
effectively, as particle data.  Thus, 3D indexing of grid data from
these datasets is not possible.

Grid Data Containers
--------------------

Data containers that return field data as multidimensional arrays
can be saved so as to preserve this type of access.  This includes
covering grids, arbitrary grids, and fixed resolution buffers.
Saving data from these containers works just as with geometric data
containers.  Field data can be accessed through geometric data
containers.

.. code-block:: python

   cg = ds.covering_grid(level=0, left_edge=[0.25]*3, dims=[16]*3)
   fn = cg.save_as_dataset(fields=["density", "particle_mass"])

   cg_ds = yt.load(fn)
   ad = cg_ds.all_data()
   print ad["grid", "density"]

Multidimensional indexing of field data is also available through
the ``data`` attribute.

.. code-block:: python

   print cg_ds.data["grid", "density"]

Fixed resolution buffers work just the same.

.. code-block:: python

   my_proj = ds.proj("density", "x", weight_field="density")
   frb = my_proj.to_frb(1.0, (800, 800))
   fn = frb.save_as_dataset(fields=["density"])
   frb_ds = yt.load(fn)
   print frb_ds.data["density"]


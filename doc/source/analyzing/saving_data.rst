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

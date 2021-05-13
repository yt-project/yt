.. _saving_data:

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
used to remake plots.  For information on this, see :ref:`remaking-plots`.

.. _saving-data-containers:

Geometric Data Containers
-------------------------

Data from geometric data containers can be saved with the
:func:`~yt.data_objects.data_containers.save_as_dataset`` function.

.. notebook-cell::

   import yt
   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")

   sphere = ds.sphere([0.5]*3, (10, "Mpc"))
   fn = sphere.save_as_dataset(fields=["density", "particle_mass"])
   print (fn)

This function will return the name of the file to which the dataset
was saved.  The filename will be a combination of the name of the
original dataset and the type of data container.  Optionally, a
specific filename can be given with the ``filename`` keyword.  If no
fields are given, the fields that have previously been queried will
be saved.

The newly created dataset can be loaded like all other supported
data through ``yt.load``.  Once loaded, field data can be accessed
through the traditional data containers or through the ``data``
attribute, which will be a data container configured like the
original data container used to make the dataset.  Grid data is
accessed by the ``grid`` data type and particle data is accessed
with the original particle type.  As with the original dataset, grid
positions and cell sizes are accessible with, for example,
("grid", "x") and ("grid", "dx").  Particle positions are
accessible as (<particle_type>, "particle_position_x").  All original
simulation parameters are accessible in the ``parameters``
dictionary, normally associated with all datasets.

.. code-block:: python

   sphere_ds = yt.load("DD0046_sphere.h5")

   # use the original data container
   print(sphere_ds.data["grid", "density"])

   # create a new data container
   ad = sphere_ds.all_data()

   # grid data
   print(ad["grid", "density"])
   print(ad["grid", "x"])
   print(ad["grid", "dx"])

   # particle data
   print(ad["all", "particle_mass"])
   print(ad["all", "particle_position_x"])

Note that because field data queried from geometric containers is
returned as unordered 1D arrays, data container datasets are treated,
effectively, as particle data.  Thus, 3D indexing of grid data from
these datasets is not possible.

.. _saving-grid-data-containers:

Grid Data Containers
--------------------

Data containers that return field data as multidimensional arrays
can be saved so as to preserve this type of access.  This includes
covering grids, arbitrary grids, and fixed resolution buffers.
Saving data from these containers works just as with geometric data
containers.  Field data can be accessed through geometric data
containers.

.. code-block:: python

   cg = ds.covering_grid(level=0, left_edge=[0.25] * 3, dims=[16] * 3)
   fn = cg.save_as_dataset(fields=["density", "particle_mass"])

   cg_ds = yt.load(fn)
   ad = cg_ds.all_data()
   print(ad["grid", "density"])

Multidimensional indexing of field data is also available through
the ``data`` attribute.

.. code-block:: python

   print(cg_ds.data["grid", "density"])

Fixed resolution buffers work just the same.

.. code-block:: python

   my_proj = ds.proj("density", "x", weight_field="density")
   frb = my_proj.to_frb(1.0, (800, 800))
   fn = frb.save_as_dataset(fields=["density"])
   frb_ds = yt.load(fn)
   print(frb_ds.data["density"])

.. _saving-spatial-plots:

Spatial Plots
-------------

Spatial plots, such as projections, slices, and off-axis slices
(cutting planes) can also be saved and reloaded.

.. code-block:: python

   proj = ds.proj("density", "x", weight_field="density")
   proj.save_as_dataset()

Once reloaded, they can be handed to their associated plotting
functions to make images.

.. code-block:: python

   proj_ds = yt.load("DD0046_proj.h5")
   p = yt.ProjectionPlot(proj_ds, "x", "density", weight_field="density")
   p.save()

.. _saving-profile-data:

Profiles
--------

Profiles created with :func:`~yt.data_objects.profiles.create_profile`,
:class:`~yt.visualization.profile_plotter.ProfilePlot`, and
:class:`~yt.visualization.profile_plotter.PhasePlot` can be saved with
the :func:`~yt.data_objects.profiles.save_as_dataset` function, which
works just as above.  Profile datasets are a type of non-spatial grid
datasets.  Geometric selection is not possible, but data can be
accessed through the ``.data`` attribute.

.. notebook-cell::

   import yt
   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   ad = ds.all_data()

   profile_2d = yt.create_profile(ad, ["density", "temperature"],
                                  "mass", weight_field=None,
                                  n_bins=(128, 128))
   profile_2d.save_as_dataset()

   prof_2d_ds = yt.load("DD0046_Profile2D.h5")
   print (prof_2d_ds.data["mass"])

The x, y (if at least 2D), and z (if 3D) bin fields can be accessed as 1D
arrays with "x", "y", and "z".

.. code-block:: python

   print(prof_2d_ds.data["x"])

The bin fields can also be returned with the same shape as the profile
data by accessing them with their original names.  This allows for
boolean masking of profile data using the bin fields.

.. code-block:: python

   # density is the x bin field
   print(prof_2d_ds.data["density"])

For 1, 2, and 3D profile datasets, a fake profile object will be
constructed by accessing the ".profile" attribute.  This is used
primarily in the case of 1 and 2D profiles to create figures using
:class:`~yt.visualization.profile_plotter.ProfilePlot` and
:class:`~yt.visualization.profile_plotter.PhasePlot`.

.. code-block:: python

   p = yt.PhasePlot(prof_2d_ds.data, "density", "temperature", "mass", weight_field=None)
   p.save()

.. _saving-array-data:

Generic Array Data
------------------

Generic arrays can be saved and reloaded as non-spatial data using
the :func:`~yt.frontends.ytdata.utilities.save_as_dataset` function,
also available as ``yt.save_as_dataset``.  As with profiles, geometric
selection is not possible, but the data can be accessed through the
``.data`` attribute.

.. notebook-cell::

   import yt
   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")

   region = ds.box([0.25]*3, [0.75]*3)
   sphere = ds.sphere(ds.domain_center, (10, "Mpc"))
   my_data = {}
   my_data["region_density"] = region["density"]
   my_data["sphere_density"] = sphere["density"]
   yt.save_as_dataset(ds, "test_data.h5", my_data)

   array_ds = yt.load("test_data.h5")
   print (array_ds.data["region_density"])
   print (array_ds.data["sphere_density"])

Array data can be saved with or without a dataset loaded.  If no
dataset has been loaded, as fake dataset can be provided as a
dictionary.

.. notebook-cell::

   import numpy as np
   import yt

   my_data = {"density": yt.YTArray(np.random.random(10), "g/cm**3"),
              "temperature": yt.YTArray(np.random.random(10), "K")}
   fake_ds = {"current_time": yt.YTQuantity(10, "Myr")}
   yt.save_as_dataset(fake_ds, "random_data.h5", my_data)

   new_ds = yt.load("random_data.h5")
   print (new_ds.data["density"])

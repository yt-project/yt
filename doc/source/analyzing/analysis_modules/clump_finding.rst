.. _clump_finding:

Clump Finding
=============

The clump finder uses a contouring algorithm to identified topologically 
disconnected structures within a dataset.  This works by first creating a 
single contour over the full range of the contouring field, then continually 
increasing the lower value of the contour until it reaches the maximum value 
of the field.  As disconnected structures are identified as separate contoures, 
the routine continues recursively through each object, creating a hierarchy of 
clumps.  Individual clumps can be kept or removed from the hierarchy based on 
the result of user-specified functions, such as checking for gravitational 
boundedness.  A sample recipe can be found in :ref:`cookbook-find_clumps`.

The clump finder requires a data object (see :ref:`data-objects`) and a field 
over which the contouring is to be performed.

.. code:: python

   import yt
   from yt.analysis_modules.level_sets.api import *

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   data_source = ds.disk([0.5, 0.5, 0.5], [0., 0., 1.],
                         (8, 'kpc'), (1, 'kpc'))

   master_clump = Clump(data_source, ("gas", "density"))

At this point, every isolated contour will be considered a clump, 
whether this is physical or not.  Validator functions can be added to 
determine if an individual contour should be considered a real clump.  
These functions are specified with the ``Clump.add_validator`` function.  
Current, two validators exist: a minimum number of cells and gravitational 
boundedness.

.. code:: python

   master_clump.add_validator("min_cells", 20)

   master_clump.add_validator("gravitationally_bound", use_particles=False)

As many validators as desired can be added, and a clump is only kept if all 
return True.  If not, a clump is remerged into its parent.  Custom validators 
can easily be added.  A validator function must only accept a ``Clump`` object 
and either return True or False.

.. code:: python

   def _minimum_gas_mass(clump, min_mass):
       return (clump["gas", "cell_mass"].sum() >= min_mass)
   add_validator("minimum_gas_mass", _minimum_gas_mass)

The ``add_validator`` function adds the validator to a registry that can 
be accessed by the clump finder.  Then, the validator can be added to the 
clump finding just like the others.

.. code:: python

   master_clump.add_validator("minimum_gas_mass", ds.quan(1.0, "Msun"))

The clump finding algorithm accepts the ``Clump`` object, the initial minimum 
and maximum of the contouring field, and the step size.  The lower value of the 
contour finder will be continually multiplied by the step size.

.. code:: python

   c_min = data_source["gas", "density"].min()
   c_max = data_source["gas", "density"].max()
   step = 2.0
   find_clumps(master_clump, c_min, c_max, step)

After the clump finding has finished, the master clump will represent the top 
of a hierarchy of clumps.  The ``children`` attribute within a ``Clump`` object 
contains a list of all sub-clumps.  Each sub-clump is also a ``Clump`` object 
with its own ``children`` attribute, and so on.

A number of helper routines exist for examining the clump hierarchy.

.. code:: python

   # Write a text file of the full hierarchy.
   write_clump_index(master_clump, 0, "%s_clump_hierarchy.txt" % ds)

   # Write a text file of only the leaf nodes.
   write_clumps(master_clump,0, "%s_clumps.txt" % ds)

   # Get a list of just the leaf nodes.
   leaf_clumps = get_lowest_clumps(master_clump)

``Clump`` objects can be used like all other data containers.

.. code:: python

   print leaf_clumps[0]["gas", "density"]
   print leaf_clumps[0].quantities.total_mass()

The writing functions will write out a series or properties about each 
clump by default.  Additional properties can be appended with the 
``Clump.add_info_item`` function.

.. code:: python

   master_clump.add_info_item("total_cells")

Just like the validators, custom info items can be added by defining functions 
that minimally accept a ``Clump`` object and return a string to be printed.

.. code:: python

   def _mass_weighted_jeans_mass(clump):
       jeans_mass = clump.data.quantities.weighted_average_quantity(
           "jeans_mass", ("gas", "cell_mass")).in_units("Msun")
       return "Jeans Mass (mass-weighted): %.6e Msolar." % jeans_mass
   add_clump_info("mass_weighted_jeans_mass", _mass_weighted_jeans_mass)

Then, add it to the list:

.. code:: python

   master_clump.add_info_item("mass_weighted_jeans_mass")

By default, the following info items are activated: **total_cells**, 
**cell_mass**, **mass_weighted_jeans_mass**, **volume_weighted_jeans_mass**, 
**max_grid_level**, **min_number_density**, **max_number_density**, and 
**distance_to_main_clump**.

Clumps can be visualized using the ``annotate_clumps`` callback.

.. code:: python

   prj = yt.ProjectionPlot(ds, 2, ("gas", "density"), 
                           center='c', width=(20,'kpc'))
   prj.annotate_clumps(leaf_clumps)
   prj.save('clumps')

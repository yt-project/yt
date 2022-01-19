.. _clump_finding:

Clump Finding
=============

The clump finder uses a contouring algorithm to identified topologically
disconnected structures within a dataset.  This works by first creating a
single contour over the full range of the contouring field, then continually
increasing the lower value of the contour until it reaches the maximum value
of the field.  As disconnected structures are identified as separate contours,
the routine continues recursively through each object, creating a hierarchy of
clumps.  Individual clumps can be kept or removed from the hierarchy based on
the result of user-specified functions, such as checking for gravitational
boundedness.  A sample recipe can be found in :ref:`cookbook-find_clumps`.

Setting up the Clump Finder
---------------------------

The clump finder requires a data object (see :ref:`data-objects`) and a field
over which the contouring is to be performed.  The data object is then used
to create the initial
:class:`~yt.data_objects.level_sets.clump_handling.Clump` object that
acts as the base for clump finding.

.. code:: python

   import yt
   from yt.data_objects.level_sets.api import *

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   data_source = ds.disk([0.5, 0.5, 0.5], [0.0, 0.0, 1.0], (8, "kpc"), (1, "kpc"))

   master_clump = Clump(data_source, ("gas", "density"))

Clump Validators
----------------

At this point, every isolated contour will be considered a clump,
whether this is physical or not.  Validator functions can be added to
determine if an individual contour should be considered a real clump.
These functions are specified with the
:func:`~yt.data_objects.level_sets.clump_handling.Clump.add_validator`
function.  Current, two validators exist: a minimum number of cells and gravitational
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
       return clump["gas", "mass"].sum() >= min_mass


   add_validator("minimum_gas_mass", _minimum_gas_mass)

The :func:`~yt.data_objects.level_sets.clump_validators.add_validator`
function adds the validator to a registry that can
be accessed by the clump finder.  Then, the validator can be added to the
clump finding just like the others.

.. code:: python

   master_clump.add_validator("minimum_gas_mass", ds.quan(1.0, "Msun"))

Running the Clump Finder
------------------------

Clump finding then proceeds by calling the
:func:`~yt.data_objects.level_sets.clump_handling.find_clumps` function.
This function accepts the
:class:`~yt.data_objects.level_sets.clump_handling.Clump` object, the initial
minimum and maximum of the contouring field, and the step size.  The lower value
of the contour finder will be continually multiplied by the step size.

.. code:: python

   c_min = data_source["gas", "density"].min()
   c_max = data_source["gas", "density"].max()
   step = 2.0
   find_clumps(master_clump, c_min, c_max, step)

Calculating Clump Quantities
----------------------------

By default, a number of quantities will be calculated for each clump when the
clump finding process has finished.  The default quantities are: ``total_cells``,
``mass``, ``mass_weighted_jeans_mass``, ``volume_weighted_jeans_mass``,
``max_grid_level``, ``min_number_density``, and ``max_number_density``.
Additional items can be added with the
:func:`~yt.data_objects.level_sets.clump_handling.Clump.add_info_item`
function.

.. code:: python

   master_clump.add_info_item("total_cells")

Just like the validators, custom info items can be added by defining functions
that minimally accept a
:class:`~yt.data_objects.level_sets.clump_handling.Clump` object and return
a format string to be printed and the value.  These are then added to the list
of available info items by calling
:func:`~yt.data_objects.level_sets.clump_info_items.add_clump_info`:

.. code:: python

   def _mass_weighted_jeans_mass(clump):
       jeans_mass = clump.data.quantities.weighted_average_quantity(
           "jeans_mass", ("gas", "mass")
       ).in_units("Msun")
       return "Jeans Mass (mass-weighted): %.6e Msolar." % jeans_mass


   add_clump_info("mass_weighted_jeans_mass", _mass_weighted_jeans_mass)

Then, add it to the list:

.. code:: python

   master_clump.add_info_item("mass_weighted_jeans_mass")

Once you have run the clump finder, you should be able to access the data for
the info item you have defined via the ``info`` attribute of a ``Clump`` object:

.. code:: python

   clump = leaf_clumps[0]
   print(clump.info["mass_weighted_jeans_mass"])

Besides the quantities calculated by default, the following are available:
``center_of_mass`` and ``distance_to_main_clump``.

Working with Clumps
-------------------

After the clump finding has finished, the master clump will represent the top
of a hierarchy of clumps.  The ``children`` attribute within a
:class:`~yt.data_objects.level_sets.clump_handling.Clump` object
contains a list of all sub-clumps.  Each sub-clump is also a
:class:`~yt.data_objects.level_sets.clump_handling.Clump` object
with its own ``children`` attribute, and so on.

.. code:: python

   print(master_clump["gas", "density"])
   print(master_clump.children)
   print(master_clump.children[0]["gas", "density"])

The entire clump tree can traversed with a loop syntax:

.. code:: python

   for clump in master_clump:
       print(clump.clump_id)

The ``leaves`` attribute of a ``Clump`` object will return a list of the
individual clumps that have no children of their own (the leaf clumps).

.. code:: python

   # Get a list of just the leaf nodes.
   leaf_clumps = master_clump.leaves

   print(leaf_clumps[0]["gas", "density"])
   print(leaf_clumps[0]["all", "particle_mass"])
   print(leaf_clumps[0].quantities.total_mass())

Visualizing Clumps
------------------

Clumps can be visualized using the ``annotate_clumps`` callback.

.. code:: python

   prj = yt.ProjectionPlot(ds, 2, ("gas", "density"), center="c", width=(20, "kpc"))
   prj.annotate_clumps(leaf_clumps)
   prj.save("clumps")

Saving and Reloading Clump Data
-------------------------------

The clump tree can be saved as a reloadable dataset with the
:func:`~yt.data_objects.level_sets.clump_handling.Clump.save_as_dataset`
function.  This will save all info items that have been calculated as well as
any field values specified with the *fields* keyword.  This function
can be called for any clump in the tree, saving that clump and all those
below it.

.. code:: python

   fn = master_clump.save_as_dataset(fields=["density", "particle_mass"])

The clump tree can then be reloaded as a regular dataset.  The ``tree`` attribute
associated with the dataset provides access to the clump tree.  The tree can be
iterated over in the same fashion as the original tree.

.. code:: python

   ds_clumps = yt.load(fn)
   for clump in ds_clumps.tree:
       print(clump.clump_id)

The ``leaves`` attribute returns a list of all leaf clumps.

.. code:: python

   print(ds_clumps.leaves)

Info items for each clump can be accessed with the ``"clump"`` field type.  Gas
or grid fields should be accessed using the ``"grid"`` field type and particle
fields should be access using the specific particle type.

.. code:: python

   my_clump = ds_clumps.leaves[0]
   print(my_clumps["clump", "mass"])
   print(my_clumps["grid", "density"])
   print(my_clumps["all", "particle_mass"])

.. _halo_catalog:

Halo Finding and Analysis
=========================

In yt-3.x, halo finding and analysis are combined into a single
framework called the
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`.
This framework is substantially different from the halo analysis
machinery available in yt-2.x and is entirely backward incompatible.
For a direct translation of various halo analysis tasks using yt-2.x
to yt-3.x, see :ref:`halo-transition`.

.. _halo_catalog_finding:

Halo Finding
------------

If you already have a halo catalog, either produced by one of the methods
below or in a format described in :ref:`halo-catalog-data`, and want to
perform further analysis, skip to :ref:`halo_catalog_analysis`.

Three halo finding methods exist within yt.  These are:

* :ref:`fof_finding`: a basic friend-of-friends algorithm (e.g. `Efstathiou et al. (1985)
  <http://adsabs.harvard.edu/abs/1985ApJS...57..241E>`_)
* :ref:`hop_finding`: `Eisenstein and Hut (1998)
  <http://adsabs.harvard.edu/abs/1998ApJ...498..137E>`_.
* :ref:`rockstar_finding`: a 6D phase-space halo finder developed by Peter Behroozi that
  scales well and does substructure finding (`Behroozi et al.
  2011 <http://adsabs.harvard.edu/abs/2011arXiv1110.4372B>`_)

Halo finding is performed through the creation of a
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`
object.  The dataset on which halo finding is to be performed should
be loaded and given to the
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`
along with the ``finder_method`` keyword to specify the method to be
used.

.. code-block:: python

   import yt
   from yt.analysis_modules.halo_analysis.api import HaloCatalog

   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, finder_method='hop')
   hc.create()

The ``finder_method`` options should be given as "fof", "hop", or
"rockstar".  Each of these methods has their own set of keyword
arguments to control functionality.  These can specified in the form
of a dictionary using the ``finder_kwargs`` keyword.

.. code-block:: python

   import yt
   from yt.analysis_modules.halo_analysis.api import HaloCatalog

   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, finder_method='fof',
                    finder_kwargs={"ptype": "stars",
                                   "padding": 0.02})
   hc.create()

For a full list of keywords for each halo finder, see
:class:`~yt.analysis_modules.halo_finding.halo_objects.FOFHaloFinder`,
:class:`~yt.analysis_modules.halo_finding.halo_objects.HOPHaloFinder`,
and
:class:`~yt.analysis_modules.halo_finding.rockstar.rockstar.RockstarHaloFinder`.

.. _fof_finding:

FOF
^^^

This is a basic friends-of-friends algorithm.  See
`Efstathiou et al. (1985)
<http://adsabs.harvard.edu/abs/1985ApJS...57..241E>`_ for more
details as well as
:class:`~yt.analysis_modules.halo_finding.halo_objects.FOFHaloFinder`.

.. _hop_finding:

HOP
^^^

The version of HOP used in yt is an upgraded version of the
`publicly available HOP code
<http://cmb.as.arizona.edu/~eisenste/hop/hop.html>`_. Support
for 64-bit floats and integers has been added, as well as
parallel analysis through spatial decomposition. HOP builds
groups in this fashion:

#. Estimates the local density at each particle using a
   smoothing kernel.

#. Builds chains of linked particles by 'hopping' from one
   particle to its densest neighbor. A particle which is
   its own densest neighbor is the end of the chain.

#. All chains that share the same densest particle are
   grouped together.

#. Groups are included, linked together, or discarded
   depending on the user-supplied over density
   threshold parameter. The default is 160.0.

See the `HOP method paper
<http://adsabs.harvard.edu/abs/1998ApJ...498..137E>`_ for
full details as well as
:class:`~yt.analysis_modules.halo_finding.halo_objects.HOPHaloFinder`.

.. _rockstar_finding:

Rockstar
^^^^^^^^

Rockstar uses an adaptive hierarchical refinement of friends-of-friends
groups in six phase-space dimensions and one time dimension, which
allows for robust (grid-independent, shape-independent, and noise-
resilient) tracking of substructure. The code is prepackaged with yt,
but also `separately available <https://bitbucket.org/gfcstanford/rockstar>`_. The lead
developer is Peter Behroozi, and the methods are described in
`Behroozi et al. 2011 <http://adsabs.harvard.edu/abs/2011arXiv1110.4372B>`_.
In order to run the Rockstar halo finder in yt, make sure you've
:ref:`installed it so that it can integrate with yt <rockstar-installation>`.

At the moment, Rockstar does not support multiple particle masses,
instead using a fixed particle mass. This will not affect most dark matter
simulations, but does make it less useful for finding halos from the stellar
mass. In simulations where the highest-resolution particles all have the
same mass (ie: zoom-in grid based simulations), one can set up a particle
filter to select the lowest mass particles and perform the halo finding
only on those.  See the this cookbook recipe for an example:
:ref:`cookbook-rockstar-nested-grid`.

To run the Rockstar Halo finding, you must launch python with MPI and
parallelization enabled. While Rockstar itself does not require MPI to run,
the MPI libraries allow yt to distribute particle information across multiple
nodes.

.. warning:: At the moment, running Rockstar inside of yt on multiple compute nodes
   connected by an Infiniband network can be problematic. Therefore, for now
   we recommend forcing the use of the non-Infiniband network (e.g. Ethernet)
   using this flag: ``--mca btl ^openib``.
   For example, here is how Rockstar might be called using 24 cores:
   ``mpirun -n 24 --mca btl ^openib python ./run_rockstar.py --parallel``.

The script above configures the Halo finder, launches a server process which
disseminates run information and coordinates writer-reader processes.
Afterwards, it launches reader and writer tasks, filling the available MPI
slots, which alternately read particle information and analyze for halo
content.

The RockstarHaloFinder class has these options that can be supplied to the
halo catalog through the ``finder_kwargs`` argument:

* ``dm_type``, the index of the dark matter particle. Default is 1.
* ``outbase``, This is where the out*list files that Rockstar makes should be
  placed. Default is 'rockstar_halos'.
* ``num_readers``, the number of reader tasks (which are idle most of the
  time.) Default is 1.
* ``num_writers``, the number of writer tasks (which are fed particles and
  do most of the analysis). Default is MPI_TASKS-num_readers-1.
  If left undefined, the above options are automatically
  configured from the number of available MPI tasks.
* ``force_res``, the resolution that Rockstar uses for various calculations
  and smoothing lengths. This is in units of Mpc/h.
  If no value is provided, this parameter is automatically set to
  the width of the smallest grid element in the simulation from the
  last data snapshot (i.e. the one where time has evolved the
  longest) in the time series:
  ``ds_last.index.get_smallest_dx() * ds_last['Mpch']``.
* ``total_particles``, if supplied, this is a pre-calculated
  total number of dark matter
  particles present in the simulation. For example, this is useful
  when analyzing a series of snapshots where the number of dark
  matter particles should not change and this will save some disk
  access time. If left unspecified, it will
  be calculated automatically. Default: ``None``.
* ``dm_only``, if set to ``True``, it will be assumed that there are
  only dark matter particles present in the simulation.
  This option does not modify the halos found by Rockstar, however
  this option can save disk access time if there are no star particles
  (or other non-dark matter particles) in the simulation. Default: ``False``.

Rockstar dumps halo information in a series of text (halo*list and
out*list) and binary (halo*bin) files inside the ``outbase`` directory.
We use the halo list classes to recover the information.

Inside the ``outbase`` directory there is a text file named ``datasets.txt``
that records the connection between ds names and the Rockstar file names.

.. _rockstar-installation:

Installing Rockstar
"""""""""""""""""""

Because of changes in the Rockstar API over time, yt only currently works with
a slightly older version of Rockstar.  This version of Rockstar has been
slightly patched and modified to run as a library inside of yt. By default it
is not installed with yt, but installation is very easy.  The
:ref:`install-script` used to install yt from source has a line:
``INST_ROCKSTAR=0`` that must be changed to ``INST_ROCKSTAR=1``.  You can
rerun this installer script over the top of an existing installation, and
it will only install components missing from the existing installation.
You can do this as follows.  Put your freshly modified install_script in
the parent directory of the yt installation directory (e.g. the parent of
``$YT_DEST``, ``yt-x86_64``, ``yt-i386``, etc.), and rerun the installer:

.. code-block:: bash

    cd $YT_DEST
    cd ..
    vi install_script.sh  // or your favorite editor to change INST_ROCKSTAR=1
    bash < install_script.sh

This will download Rockstar and install it as a library in yt.

.. _halo_catalog_analysis:

Extra Halo Analysis
-------------------

As a reminder, all halo catalogs created by the methods outlined in
:ref:`halo_catalog_finding` as well as those in the formats discussed in
:ref:`halo-catalog-data` can be loaded in to yt as first-class datasets.
Once a halo catalog has been created, further analysis can be performed
by providing both the halo catalog and the original simulation dataset to
the
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`.

.. code-block:: python

   halos_ds = yt.load('rockstar_halos/halos_0.0.bin')
   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, halos_ds=halos_ds)

A data object can also be supplied via the keyword ``data_source``,
associated with either dataset, to control the spatial region in
which halo analysis will be performed.

The :class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`
allows the user to create a pipeline of analysis actions that will be
performed on all halos in the existing catalog.  The analysis can be
performed in parallel with separate processors or groups of processors
being allocated to perform the entire pipeline on individual halos.
The pipeline is setup by adding actions to the
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`.
Each action is represented by a callback function that will be run on
each halo.  There are four types of actions:

* :ref:`halo_catalog_filters`
* :ref:`halo_catalog_quantities`
* :ref:`halo_catalog_callbacks`
* :ref:`halo_catalog_recipes`

A list of all available filters, quantities, and callbacks can be found in
:ref:`halo_analysis_ref`.
All interaction with this analysis can be performed by importing from
halo_analysis.

.. _halo_catalog_filters:

Filters
^^^^^^^

A filter is a function that returns True or False. If the return value
is True, any further queued analysis will proceed and the halo in
question will be added to the final catalog. If the return value False,
further analysis will not be performed and the halo will not be included
in the final catalog.

An example of adding a filter:

.. code-block:: python

   hc.add_filter('quantity_value', 'particle_mass', '>', 1E13, 'Msun')

Currently quantity_value is the only available filter, but more can be
added by the user by defining a function that accepts a halo object as
the first argument and then adding it as an available filter. If you
think that your filter may be of use to the general community, you can
add it to ``yt/analysis_modules/halo_analysis/halo_filters.py`` and issue a
pull request.

An example of defining your own filter:

.. code-block:: python

   def my_filter_function(halo):

       # Define condition for filter
       filter_value = True

       # Return a boolean value
       return filter_value

   # Add your filter to the filter registry
   add_filter("my_filter", my_filter_function)

   # ... Later on in your script
   hc.add_filter("my_filter")

.. _halo_catalog_quantities:

Quantities
^^^^^^^^^^

A quantity is a call back that returns a value or values. The return values
are stored within the halo object in a dictionary called "quantities." At
the end of the analysis, all of these quantities will be written to disk as
the final form of the generated halo catalog.

Quantities may be available in the initial fields found in the halo catalog,
or calculated from a function after supplying a definition. An example
definition of center of mass is shown below. Currently available quantities
are center_of_mass and bulk_velocity. Their definitions are available in
``yt/analysis_modules/halo_analysis/halo_quantities.py``. If you think that
your quantity may be of use to the general community, add it to
``halo_quantities.py`` and issue a pull request.  Default halo quantities are:

* ``particle_identifier`` -- Halo ID (e.g. 0 to N)
* ``particle_mass`` -- Mass of halo
* ``particle_position_x`` -- Location of halo
* ``particle_position_y`` -- Location of halo
* ``particle_position_z`` -- Location of halo
* ``virial_radius`` -- Virial radius of halo

An example of adding a quantity:

.. code-block:: python

   hc.add_quantity('center_of_mass')

An example of defining your own quantity:

.. code-block:: python

   def my_quantity_function(halo):
       # Define quantity to return
       quantity = 5

       return quantity

   # Add your filter to the filter registry
   add_quantity('my_quantity', my_quantity_function)


   # ... Later on in your script
   hc.add_quantity("my_quantity")

This quantity will then be accessible for functions called later via the
*quantities* dictionary that is associated with the halo object.

.. code-block:: python

   def my_new_function(halo):
       print(halo.quantities["my_quantity"])
   add_callback("print_quantity", my_new_function)

   # ... Anywhere after "my_quantity" has been called
   hc.add_callback("print_quantity")

.. _halo_catalog_callbacks:

Callbacks
^^^^^^^^^

A callback is actually the super class for quantities and filters and
is a general purpose function that does something, anything, to a Halo
object. This can include hanging new attributes off the Halo object,
performing analysis and writing to disk, etc. A callback does not return
anything.

An example of using a pre-defined callback where we create a sphere for
each halo with a radius that is twice the saved ``radius``.

.. code-block:: python

   hc.add_callback("sphere", factor=2.0)

Currently available callbacks are located in
``yt/analysis_modules/halo_analysis/halo_callbacks.py``.  New callbacks may
be added by using the syntax shown below. If you think that your
callback may be of use to the general community, add it to
halo_callbacks.py and issue a pull request.

An example of defining your own callback:

.. code-block:: python

   def my_callback_function(halo):
       # Perform some callback actions here
       x = 2
       halo.x_val = x

   # Add the callback to the callback registry
   add_callback('my_callback', my_callback_function)


   # ...  Later on in your script
   hc.add_callback("my_callback")

.. _halo_catalog_recipes:

Recipes
^^^^^^^

Recipes allow you to create analysis tasks that consist of a series of
callbacks, quantities, and filters that are run in succession.  An example
of this is
:func:`~yt.analysis_modules.halo_analysis.halo_recipes.calculate_virial_quantities`,
which calculates virial quantities by first creating a sphere container,
performing 1D radial profiles, and then interpolating to get values at a
specified threshold overdensity.  All of these operations are separate
callbacks, but the recipes allow you to add them to your analysis pipeline
with one call.  For example,

.. code-block:: python

   hc.add_recipe("calculate_virial_quantities", ["radius", "matter_mass"])

The available recipes are located in
``yt/analysis_modules/halo_analysis/halo_recipes.py``.  New recipes can be
created in the following manner:

.. code-block:: python

   def my_recipe(halo_catalog, fields, weight_field=None):
       # create a sphere
       halo_catalog.add_callback("sphere")
       # make profiles
       halo_catalog.add_callback("profile", ["radius"], fields,
                                 weight_field=weight_field)
       # save the profile data
       halo_catalog.add_callback("save_profiles", output_dir="profiles")

   # add recipe to the registry of recipes
   add_recipe("profile_and_save", my_recipe)


   # ...  Later on in your script
   hc.add_recipe("profile_and_save", ["density", "temperature"],
                 weight_field="cell_mass")

Note, that unlike callback, filter, and quantity functions that take a ``Halo``
object as the first argument, recipe functions should take a ``HaloCatalog``
object as the first argument.

Running the Pipeline
--------------------

After all callbacks, quantities, and filters have been added, the
analysis begins with a call to HaloCatalog.create.

.. code-block:: python

   hc.create()

The save_halos keyword determines whether the actual Halo objects
are saved after analysis on them has completed or whether just the
contents of their quantities dicts will be retained for creating the
final catalog. The looping over halos uses a call to parallel_objects
allowing the user to control how many processors work on each halo.
The final catalog is written to disk in the output directory given
when the
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`
object was created.

All callbacks, quantities, and filters are stored in an actions list,
meaning that they are executed in the same order in which they were added.
This enables the use of simple, reusable, single action callbacks that
depend on each other. This also prevents unnecessary computation by allowing
the user to add filters at multiple stages to skip remaining analysis if it
is not warranted.

Saving and Reloading Halo Catalogs
----------------------------------

A :class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`
saved to disk can be reloaded as a yt dataset with the
standard call to ``yt.load``.  See :ref:`halocatalog` for a demonstration
of loading and working only with the catalog.
Any side data, such as profiles, can be reloaded
with a ``load_profiles`` callback and a call to
:func:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog.load`.

.. code-block:: python

   hds = yt.load(path+"halo_catalogs/catalog_0046/catalog_0046.0.h5")
   hc = HaloCatalog(halos_ds=hds,
                    output_dir="halo_catalogs/catalog_0046")
   hc.add_callback("load_profiles", output_dir="profiles",
                   filename="virial_profiles")
   hc.load()

Halo Catalog in Action
----------------------

For a full example of how to use these methods together see
:ref:`halo-analysis-example`.

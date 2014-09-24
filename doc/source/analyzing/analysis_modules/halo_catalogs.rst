.. _halo_catalog:

Halo Catalogs
=============

Creating Halo Catalogs
----------------------

In yt 3.0, operations relating to the analysis of halos (halo finding,
merger tree creation, and individual halo analysis) are all brought 
together into a single framework. This framework is substantially
different from the halo analysis machinery available in yt-2.x and is 
entirely backward incompatible.  
For a direct translation of various halo analysis tasks using yt-2.x
to yt-3.0 please see :ref:`halo-transition`.

A catalog of halos can be created from any initial dataset given to halo 
catalog through data_ds. These halos can be found using friends-of-friends,
HOP, and Rockstar. The finder_method keyword dictates which halo finder to
use. The available arguments are :ref:`fof`, :ref:`hop`, and :ref:`rockstar`. 
For more details on the relative differences between these halo finders see 
:ref:`halo_finding`.

The class which holds all of the halo information is the 
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`.

.. code-block:: python

   import yt
   from yt.analysis_modules.halo_analysis.api import HaloCatalog

   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, finder_method='hop')

A halo catalog may also be created from already run rockstar outputs. 
This method is not implemented for previously run friends-of-friends or 
HOP finders. Even though rockstar creates one file per processor, 
specifying any one file allows the full catalog to be loaded. Here we 
only specify the file output by the processor with ID 0. Note that the 
argument for supplying a rockstar output is `halos_ds`, not `data_ds`.

.. code-block:: python

   halos_ds = yt.load(path+'rockstar_halos/halos_0.0.bin')
   hc = HaloCatalog(halos_ds=halos_ds)

Although supplying only the binary output of the rockstar halo finder 
is sufficient for creating a halo catalog, it is not possible to find 
any new information about the identified halos. To associate the halos 
with the dataset from which they were found, supply arguments to both 
halos_ds and data_ds.

.. code-block:: python

   halos_ds = yt.load(path+'rockstar_halos/halos_0.0.bin')
   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, halos_ds=halos_ds)

A data object can also be supplied via the keyword ``data_source``, 
associated with either dataset, to control the spatial region in 
which halo analysis will be performed.

Analysis Using Halo Catalogs
----------------------------

Analysis is done by adding actions to the 
:class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog`.
Each action is represented by a callback function that will be run on each halo. 
There are three types of actions:

* Filters
* Quantities
* Callbacks

A list of all available filters, quantities, and callbacks can be found in 
:ref:`halo_analysis_ref`.  
All interaction with this analysis can be performed by importing from 
halo_analysis.

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

Quantities
^^^^^^^^^^

A quantity is a call back that returns a value or values. The return values 
are stored within the halo object in a dictionary called “quantities.” At 
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
       print halo.quantities["my_quantity"]
   add_callback("print_quantity", my_new_function)

   # ... Anywhere after "my_quantity" has been called
   hc.add_callback("print_quantity")

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

Running Analysis
----------------

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
depend on each other. This also prevents unecessary computation by allowing 
the user to add filters at multiple stages to skip remaining analysis if it 
is not warranted.

Saving and Reloading Halo Catalogs
----------------------------------

A :class:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog` 
saved to disk can be reloaded as a yt dataset with the 
standard call to load. Any side data, such as profiles, can be reloaded 
with a ``load_profiles`` callback and a call to 
:func:`~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog.load`.

.. code-block:: python

   hds = yt.load(path+"halo_catalogs/catalog_0046/catalog_0046.0.h5")
   hc = HaloCatalog(halos_ds=hds,
                    output_dir="halo_catalogs/catalog_0046")
   hc.add_callback("load_profiles", output_dir="profiles",
                   filename="virial_profiles")
   hc.load()

Worked Example of Halo Catalog in Action
----------------------------------------

For a full example of how to use these methods together see 
:ref:`halo-analysis-example`.

.. _halo-transition:

Getting up to Speed with Halo Analysis in yt-3.0
================================================

If you're used to halo analysis in yt-2.x, heres a guide to
how to update your analysis pipeline to take advantage of
the new halo catalog infrastructure. 

Finding Halos
-------------

Previously, halos were found using calls to ``HaloFinder``, 
``FOFHaloFinder`` and ``RockstarHaloFinder``. Now it is 
encouraged that you find the halos upon creation of the halo catalog 
by supplying a value to the ``finder_method`` keyword when calling
``HaloCatalog``. Currently, only halos found using rockstar or a 
previous instance of a halo catalog are able to be loaded 
using the ``halos_ds`` keyword.

To pass additional arguments to the halo finders 
themselves, supply a dictionary to ``finder_kwargs`` where
each key in the dictionary is a keyword of the halo finder
and the corresponding value is the value to be passed for
that keyword.

Getting Halo Information
------------------------
All quantities that used to be present in a ``halo_list`` are
still able to be found but are not necessarily included by default.
Every halo will by default have the following properties:

* particle_position_i (where i can be x,y,z)
* particle_mass
* virial_radius
* particle_identifier

If other quantities are desired, they can be included by adding
the corresponding quantity before the catalog is created. See
the full halo catalog documentation for further information about
how to add these quantities and what quantities are available.

You no longer have to iteratre over halos in the ``halo_list``.
Now a halo dataset can be treated as a regular dataset and 
all quantities are available by accessing ``all_data``.
Specifically, all quantities can be accessed as shown:

.. code-block:: python

   import yt
   from yt.analysis_modules.halo_analysis.api import HaloCatalog
   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, finder_method='hop')
   hc.create()
   ad = hc.halos_ds.all_data()
   masses = ad['particle_mass'][:]


Prefiltering Halos
------------------

Prefiltering halos before analysis takes place is now done
by adding a filter before the call to create. An example
is shown below

.. code-block:: python

   import yt
   from yt.analysis_modules.halo_analysis.api import HaloCatalog
   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, finder_method='hop')
   hc.add_filter("quantity_value", "particle_mass", ">", 1e13, "Msun")
   hc.create()

Profiling Halos
---------------

The halo profiler available in yt-2.x has been removed, and
profiling functionality is now completely contained within the
halo catalog. A complete example of how to profile halos by 
radius using the new infrastructure is given in 
:ref:`halo-analysis-example`. 

Plotting Halos
--------------

Annotating halo locations onto a slice or projection works in 
the same way as in yt-2.x, but now a halo catalog must be
passed to the annotate halo call rather than a halo list.

.. code-block:: python

   import yt
   from yt.analysis_modules.halo_analysis.api import HaloCatalog

   data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
   hc = HaloCatalog(data_ds=data_ds, finder_method='hop')
   hc.create()

   prj = yt.ProjectionPlot(data_ds, 'z', 'density')
   prj.annotate_halos(hc)
   prj.save()

Written Data
------------

Data is now written out in the form of h5 files rather than
text files. The directory they are written out to is 
controlled by the keyword ``output_dir``. Each quantity
is a field in the file.

.. _light-ray-generator:

Light Ray Generator
====================
.. sectionauthor:: Britton Smith <brittonsmith@gmail.com>

Light rays are similar to light cones (:ref:`light-cone-generator`) in how  
they stack multiple datasets together to span a redshift interval.  Unlike 
light cones, which which stack randomly oriented projections from each 
dataset to create synthetic images, light rays use thin pencil beams to 
simulate QSO sight lines.

.. image:: _images/lightray.png

A ray segment records the information of all grid cells intersected by the ray 
as well as the path length, dl, of the ray through the cell.  Column densities 
can be calculated by multiplying physical densities by the path length.

Configuring the Light Ray Generator
-----------------------------------
  
The arguments required to instantiate a ``LightRay`` object are the same as 
those required for a ``LightCone`` object: the simulation parameter file, the 
simulation type, the nearest redshift, and the furthest redshift.

.. code-block:: python

  from yt.analysis_modules.cosmological_observation.api import LightRay
  lr = LightRay("enzo_tiny_cosmology/32Mpc_32.enzo",
                'Enzo', 0.0, 0.1)

Additional keyword arguments are:

 * **use_minimum_datasets** (*bool*): If True, the minimum number of datasets 
   is used to connect the initial and final redshift.  If false, the light 
   ray solution will contain as many entries as possible within the redshift
   interval.  Default: True.

 * **deltaz_min** (*float*):  Specifies the minimum Delta-z between consecutive
   datasets in the returned list.  Default: 0.0.

 * **minimum_coherent_box_fraction** (*float*): Used with use_minimum_datasets 
   set to False, this parameter specifies the fraction of the total box size 
   to be traversed before rerandomizing the projection axis and center.  This
   was invented to allow light rays with thin slices to sample coherent large 
   scale structure, but in practice does not work so well.  Try setting this 
   parameter to 1 and see what happens.  Default: 0.0.

 * **time_data** (*bool*): Whether or not to include time outputs when gathering
   datasets for time series.  Default: True.

 * **redshift_data** (*bool*): Whether or not to include redshift outputs when 
   gathering datasets for time series.  Default: True.


Making Light Ray Data
---------------------

Once the LightRay object has been instantiated, the :meth:`make_light_ray` 
will trace out the rays in each dataset and collect information for all the 
fields requested.  The output file will be an hdf5 file containing all the 
cell field values for all the cells that were intersected by the ray.  A 
single LightRay object can be used over and over to make multiple 
randomizations, simply by changing the value of the random seed with the 
**seed** keyword.

.. code-block:: python

  lr.make_light_ray(seed=8675309,
                    fields=['temperature', 'density'],
                    get_los_velocity=True)

The keyword arguments are:

 * **seed** (*int*): Seed for the random number generator.  Default: None.

 * **fields** (*list*): A list of fields for which to get data.  Default: None.

 * **solution_filename** (*string*): Path to a text file where the 
   trajectories of each subray is written out.  Default: None.

 * **data_filename** (*string*): Path to output file for ray data.  
   Default: None.

 * **get_los_velocity** (*bool*): If True, the line of sight velocity is 
   calculated for each point in the ray.  Default: False.

 * **get_nearest_halo** (*bool*): If True, the HaloProfiler will be used to 
   calculate the distance and mass of the nearest halo for each point in the
   ray.  This option requires additional information to be included.  See 
   the cookbook for an example.  Default: False.

 * **nearest_halo_fields** (*list*): A list of fields to be calculated for the 
   halos nearest to every pixel in the ray.  Default: None.

 * **halo_list_file** (*str*): Filename containing a list of halo properties to be used 
   for getting the nearest halos to absorbers.  Default: None.

 * **halo_profiler_parameters** (*dict*): A dictionary of parameters to be 
   passed to the HaloProfiler to create the appropriate data used to get 
   properties for the nearest halos.  Default: None.

 * **njobs** (*int*): The number of parallel jobs over which the slices for the
   halo mask will be split.  Choose -1 for one processor per individual slice 
   and 1 to have all processors work together on each projection.  Default: 1

 * **dynamic** (*bool*): If True, use dynamic load balancing to create the 
   projections.  Default: False.

Getting The Nearest Galaxies
----------------------------

The light ray tool will use the HaloProfiler to calculate the distance and 
mass of the nearest halo to that pixel.  In order to do this, a dictionary 
called halo_profiler_parameters is used to pass instructions to the 
HaloProfiler.  This dictionary has three additional keywords:

 * **halo_profiler_kwargs** (*dict*): A dictionary of standard HaloProfiler 
   keyword arguments and values to be given to the HaloProfiler.

 * **halo_profiler_actions** (*list*): A list of actions to be performed by 
   the HaloProfiler.  Each item in the list should be a dictionary with the 
   following entries: "function", "args", and "kwargs", for the function to 
   be performed, the arguments supplied to that function, and the keyword 
   arguments.

 * **halo_list** (*string*): 'all' to use the full halo list, or 'filtered' 
   to use the filtered halo list created after calling make_profiles.

See the recipe in the cookbook for am example.

What Can I do with this?
------------------------

Try :ref:`absorption_spectrum`.
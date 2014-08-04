.. _light-ray-generator:

Light Ray Generator
===================

Light rays are similar to light cones (:ref:`light-cone-generator`) in how  
they stack multiple datasets together to span a redshift interval.  Unlike 
light cones, which stack randomly oriented projections from each 
dataset to create synthetic images, light rays use thin pencil beams to 
simulate QSO sight lines.  A sample script can be found in the cookbook 
under :ref:`cookbook-light_ray`.

.. image:: _images/lightray.png

A ray segment records the information of all grid cells intersected by the 
ray as well as the path length, ``dl``, of the ray through the cell.  Column 
densities can be calculated by multiplying physical densities by the path 
length.

Configuring the Light Ray Generator
-----------------------------------
  
The arguments required to instantiate a 
:class:`~yt.analysis_modules.cosmological_observation.light_ray.light_ray.LightRay` 
object are the same as 
those required for a 
:class:`~yt.analysis_modules.cosmological_observation.light_cone.light_cone.LightCone` 
object: the simulation parameter file, the 
simulation type, the nearest redshift, and the furthest redshift.

.. code-block:: python

  from yt.analysis_modules.cosmological_observation.api import LightRay
  lr = LightRay("enzo_tiny_cosmology/32Mpc_32.enzo",
                'Enzo', 0.0, 0.1)

Additional keyword arguments are:

* ``use_minimum_datasets`` (*bool*): If True, the minimum number of datasets 
  is used to connect the initial and final redshift.  If false, the light 
  ray solution will contain as many entries as possible within the redshift
  interval.  Default: True.

* ``deltaz_min`` (*float*):  Specifies the minimum Delta-z between 
  consecutive datasets in the returned list.  Default: 0.0.

* ``minimum_coherent_box_fraction`` (*float*): Used with 
  ``use_minimum_datasets`` set to False, this parameter specifies the 
  fraction of the total box size to be traversed before rerandomizing the 
  projection axis and center.  This was invented to allow light rays with 
  thin slices to sample coherent large scale structure, but in practice 
  does not work so well.  Try setting this parameter to 1 and see what 
  happens.  Default: 0.0.

* ``time_data`` (*bool*): Whether or not to include time outputs when 
  gathering datasets for time series.  Default: True.

* ``redshift_data`` (*bool*): Whether or not to include redshift outputs 
  when gathering datasets for time series.  Default: True.

Making Light Ray Data
---------------------

Once the LightRay object has been instantiated, the 
:func:`~yt.analysis_modules.cosmological_observation.light_ray.light_ray.LightRay,make_light_ray` 
function will trace out the rays in each dataset and collect information for all the 
fields requested.  The output file will be an HDF5 file containing all the 
cell field values for all the cells that were intersected by the ray.  A 
single LightRay object can be used over and over to make multiple 
randomizations, simply by changing the value of the random seed with the 
``seed`` keyword.

.. code-block:: python

  lr.make_light_ray(seed=8675309,
                    fields=['temperature', 'density'],
                    get_los_velocity=True)

The keyword arguments are:

* ``seed`` (*int*): Seed for the random number generator.  Default: None.

* ``start_position`` (*list* of floats): Used only if creating a light ray 
  from a single dataset.  The coordinates of the starting position of the 
  ray.  Default: None.

* ``end_position`` (*list* of floats): Used only if creating a light ray 
  from a single dataset.  The coordinates of the ending position of the ray.
  Default: None.

* ``trajectory`` (*list* of floats): Used only if creating a light ray 
  from a single dataset.  The (r, theta, phi) direction of the light ray.  
  Use either ``end_position`` or ``trajectory``, not both.  
  Default: None.

* ``fields`` (*list*): A list of fields for which to get data.  
  Default: None.

* ``solution_filename`` (*string*): Path to a text file where the 
  trajectories of each subray is written out.  Default: None.

* ``data_filename`` (*string*): Path to output file for ray data.  
  Default: None.

* ``get_los_velocity`` (*bool*): If True, the line of sight velocity is 
  calculated for each point in the ray.  Default: True.

* ``njobs`` (*int*): The number of parallel jobs over which the slices for 
  the halo mask will be split.  Choose -1 for one processor per individual 
  slice and 1 to have all processors work together on each projection.  
  Default: 1

* ``dynamic`` (*bool*): If True, use dynamic load balancing to create the 
  projections.  Default: False.

.. note:: As of :code:`yt-3.0`, the functionality for recording properties of the nearest halo to each element of the ray no longer exists.  This is still available in :code:`yt-2.x`.  If you would like to use this feature in :code:`yt-3.x`, help is needed to port it over.  Contact the yt-users mailing list if you are interested in doing this.

What Can I do with this?
------------------------

Try :ref:`absorption_spectrum`.

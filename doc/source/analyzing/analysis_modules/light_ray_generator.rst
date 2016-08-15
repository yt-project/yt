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

Below follows the creation of a light ray from multiple datasets stacked
together.  However, a light ray can also be made from a single dataset.
For an example of this, see :ref:`cookbook-single-dataset-light-ray`.

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
                simulation_type="Enzo",
                near_redshift=0.0, far_redshift=0.1)

Additional keyword arguments are:

* ``use_minimum_datasets`` (*bool*): If True, the minimum number of datasets
  is used to connect the initial and final redshift.  If false, the light
  ray solution will contain as many entries as possible within the redshift
  interval.  Default: True.

* ``deltaz_min`` (*float*):  Specifies the minimum Delta-z between
  consecutive datasets in the returned list.  Default: 0.0.

* ``minimum_coherent_box_fraction`` (*float*): Use to specify the minimum
  length of a ray, in terms of the size of the domain, before the trajectory
  is re-randomized.  Set to 0 to have ray trajectory randomized for every
  dataset.  Set to np.inf (infinity) to use a single trajectory for the
  entire ray.  Default: 0.0.

* ``high_res_box_size_fraction`` (*float*): For use with zoom-in simulations,
  use to specify the size of the high resolution region of the simulation in
  terms of the fraction of the total domain size.  If set, the light ray
  solution will be calculated such that rays only make use of the high
  resolution region.  Default: 1.0.

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
                    use_peculiar_velocity=True)

The keyword arguments are:

* ``seed`` (*int*): Seed for the random number generator.  Default: None.

* ``periodic`` (*bool*): If True, ray trajectories will make use of periodic
  boundaries.  If False, ray trajectories will not be periodic.  Default : True.

* ``left_edge`` (iterable of *floats* or *YTArray*): The left corner of the
  region in which rays are to be generated.  If None, the left edge will be
  that of the domain.  Default: None.

* ``right_edge`` (iterable of *floats* or *YTArray*): The right corner of
  the region in which rays are to be generated.  If None, the right edge
  will be that of the domain.  Default: None.

* ``min_level`` (*int*): The minimum refinement level of the spatial region in
  which the ray passes.  This can be used with zoom-in simulations where the
  high resolution region does not keep a constant geometry.  Default: None.

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

* ``use_peculiar_velocity`` (*bool*): If True, the doppler redshift from
  the peculiar velocity of gas along the ray is calculated and added to the
  cosmological redshift as the "effective" redshift.
  Default: True.

* ``redshift`` (*float*): Used with light rays made from single datasets to
  specify a starting redshift for the ray.  If not used, the starting
  redshift will be 0 for a non-cosmological dataset and the dataset redshift
  for a cosmological dataset.  Default: None.

* ``njobs`` (*int*): The number of parallel jobs over which the slices for
  the halo mask will be split.  Choose -1 for one processor per individual
  slice and 1 to have all processors work together on each projection.
  Default: 1

.. note:: As of :code:`yt-3.0`, the functionality for recording properties of the nearest halo to each element of the ray no longer exists.  This is still available in :code:`yt-2.x`.  If you would like to use this feature in :code:`yt-3.x`, help is needed to port it over.  Contact the yt-users mailing list if you are interested in doing this.

What Can I do with this?
------------------------

Once you have created a `LightRay`, you can use it to generate an
:ref:`absorption_spectrum`.  In addition, you can use the
:class:`~yt.visualization.plot_modifications.RayCallback` to
:ref:`annotate-ray` on your plots.

.. _light-cone-generator:

Light Cone Generator
====================
.. sectionauthor:: Britton Smith <brittonsmith@gmail.com>

Light cones are projections made by stacking multiple datasets together to 
continuously span a given redshift interval.  The width of individual 
projection slices is adjusted such that each slice has the same angular size.  
Each projection slice is randomly shifted and projected along a random axis to 
ensure that the same structures are not sampled multiple times.  Since deeper 
images sample earlier epochs of the simulation, light cones represent the 
closest thing to synthetic imaging observations.

.. image:: _images/LightCone_full_small.png
   :width: 500

A light cone projection of the thermal Sunyaev-Zeldovich Y parameter from 
z = 0 to 0.4 with a 450x450 arcminute field of view using 9 individual 
slices.  The panels shows the contributions from the 9 individual slices with 
the final light cone image shown in the bottom, right.

Configuring the Light Cone Generator
------------------------------------

A recipe for creating a simple light cone projection can be found in the 
cookbook.  The required arguments to instantiate a ``LightCone`` objects are 
the path to the simulation parameter file, the simulation type, the nearest 
redshift, and the furthest redshift of the light cone.

.. code-block:: python

  from yt.analysis_modules.api import LightCone

  lc = LightCone('enzo_tiny_cosmology/32Mpc_32.enzo',
                 'Enzo', 0., 0.1)

The additional keyword arguments are:

 * **field_of_view_in_arcminutes** (*float*): The field of view of the image 
   in units of arcminutes.  Default: 600.0.

 * **image_resolution_in_arcseconds** (*float*): The size of each image pixel 
   in units of arcseconds.  Default: 60.0.

 * **use_minimum_datasets** (*bool*):  If True, the minimum number of datasets 
   is used to connect the initial and final redshift.  If false, the light 
   cone solution will contain as many entries as possible within the redshift 
   interval.  Default: True.

 * **deltaz_min** (*float*): Specifies the minimum Delta-z between 
   consecutive datasets in the returned list.  Default: 0.0.

 * **minimum_coherent_box_fraction** (*float*): Used with use_minimum_datasets 
   set to False, this parameter specifies the fraction of the total box size 
   to be traversed before rerandomizing the projection axis and center.  This 
   was invented to allow light cones with thin slices to sample coherent large 
   scale structure, but in practice does not work so well.  Try setting this 
   parameter to 1 and see what happens.  Default: 0.0.

 * **time_data** (*bool*): Whether or not to include time outputs when 
   gathering datasets for time series.  Default: True.

 * **redshift_data** (*bool*): Whether or not to include redshift outputs when 
   gathering datasets for time series.  Default: True.

 * **set_parameters** (*dict*): Dictionary of parameters to attach to 
   ds.parameters.  Default: None.

 * **output_dir** (*string*): The directory in which images and data files
    will be written.  Default: 'LC'.

 * **output_prefix** (*string*): The prefix of all images and data files.
   Default: 'LightCone'.

Creating Light Cone Solutions
-----------------------------

A light cone solution consists of a list of datasets and the width, depth, 
center, and axis of the projection to be made for that slice.  The 
:meth:`LightCone.calculate_light_cone_solution` function is used to 
calculate the random shifting and projection axis:

.. code-block:: python

  lc.calculate_light_cone_solution(seed=123456789, filename='lightcone.dat')

The keyword argument are:

 * **seed** (*int*): the seed for the random number generator.  Any light cone 
   solution can be reproduced by giving the same random seed.  Default: None 
   (each solution will be distinct).

 * **filename** (*str*): if given, a text file detailing the solution will be 
   written out.  Default: None.

If a new solution for the same LightCone object is desired, the 
:meth:`rerandomize_light_cone_solution` method should be called in place of 
:meth:`calculate_light_cone_solution`:

.. code-block:: python

  new_seed = 987654321
  lc.rerandomize_light_cone_solution(new_seed, Recycle=True, 
                                     filename='new_lightcone.dat')

Additional keyword arguments are:

 * **recycle** (*bool*): if True, the new solution will have the same shift in 
   the line of sight as the original solution.  Since the projections of each 
   slice are serialized and stored for the entire width of the box (even if 
   the width used is left than the total box), the projection data can be 
   deserialized instead of being remade from scratch.  This can greatly speed 
   up the creation of a large number of light cone projections.  Default: True.

 * **filename** (*str*): if given, a text file detailing the solution will be 
   written out.  Default: None.

If :meth:`rerandomize_light_cone_solution` is used, the LightCone object will 
keep a copy of the original solution that can be returned to at any time by 
calling :meth:`restore_master_solution`:

.. code-block:: python

  lc.restore_master_solution()

.. note:: All light cone solutions made with the above method will still use 
   the same list of datasets.  Only the shifting and projection axis will be 
   different.

Making a Light Cone Projection
------------------------------

With the light cone solution set, projections can be made of any available 
field:

.. code-block:: python

  field = 'density'
  lc.project_light_cone(field , weight_field=None, 
                        save_stack=True, 
                        save_slice_images=True)

Additional keyword arguments:

 * **weight_field** (*str*): the weight field of the projection.  This has the 
   same meaning as in standard projections.  Default: None.

 * **apply_halo_mask** (*bool*): if True, a boolean mask is apply to the light 
   cone projection.  See below for a description of halo masks.  Default: False.

 * **node** (*str*): a prefix to be prepended to the node name under which the 
   projection data is serialized.  Default: None.

 * **save_stack** (*bool*): if True, the unflatted light cone data including 
   each individual slice is written to an hdf5 file.  Default: True.

 * **save_final_image** (*bool*): if True, save an image of the final light 
   cone projection.  Default: True.

 * **save_slice_images** (*bool*): save images for each individual projection 
   slice.  Default: False.

 * **flatten_stack** (*bool*): if True, the light cone stack is continually 
   flattened each time a slice is added in order to save memory.  This is 
   generally not necessary.  Default: False.

 * **photon_field** (*bool*): if True, the projection data for each slice is 
   decremented by 4 pi R :superscript:`2` , where R is the luminosity 
   distance between the observer and the slice redshift.  Default: False.

 * **njobs** (*int*): The number of parallel jobs over which the light cone 
   projection will be split.  Choose -1 for one processor per individual
   projection and 1 to have all processors work together on each projection.
   Default: 1.

 * **dynamic** (*bool*): If True, use dynamic load balancing to create the 
   projections.  Default: False.

Sampling Unique Light Cone Volumes
----------------------------------

When making a large number of light cones, particularly for statistical 
analysis, it is important to have a handle on the amount of sampled volume in 
common from one projection to another.  Any statistics may untrustworthy if a 
set of light cones have too much volume in common, even if they may all be 
entirely different in appearance.  LightCone objects have the ability to 
calculate the volume in common between two solutions with the same dataset 
ist.  The :meth:`find_unique_solutions` and 
:meth:`project_unique_light_cones` functions can be used to create a set of 
light cone solutions that have some maximum volume in common and create light 
cone projections for those solutions.  If specified, the code will attempt to 
use recycled solutions that can use the same serialized projection objects 
that have already been created.  This can greatly increase the speed of making 
multiple light cone projections.  See the cookbook for an example of doing this.

Making Light Cones with a Halo Mask
-----------------------------------

The situation may arise where it is necessary or desirable to know the 
location of halos within the light cone volume, and specifically their 
location in the final image.  This can be useful for developing algorithms to 
find galaxies or clusters in image data.  The light cone generator does this 
by running the HaloProfiler (see :ref:`halo_profiling`) on each of the 
datasets used in the light cone and shifting them accordingly with the light 
cone solution.  The ability also exists to create a boolean mask with the 
dimensions of the final light cone image that can be used to mask out the 
halos in the image.  It is left as an exercise to the reader to find a use for 
this functionality.  This process is somewhat complicated, but not terribly.  
See the recipe in the cookbook for an example of this functionality.

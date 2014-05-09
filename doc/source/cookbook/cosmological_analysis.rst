Cosmological Analysis
---------------------

These scripts demonstrate some basic and more advanced analysis that can be 
performed on cosmological simulations.

Plotting Halos
~~~~~~~~~~~~~~
This is a mechanism for plotting circles representing identified particle halos
on an image.

.. yt_cookbook:: halo_plotting.py

.. _cookbook-halo_finding:

Halo Profiling and Custom Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This script demonstrates the three primary uses of the halo profiler: 
1) radial profiles and filtering; 2) projections; and 3) custom halo 
analysis.

.. yt_cookbook:: halo_profiler.py

Halo Tracking Across Timesteps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This script demonstrates tracking a halo across multiple timesteps
in a TimeSeries object, as well as some handy functions for looking
at the properties of that halo over time.

.. yt_cookbook:: halo_merger_tree.py

Light Cone Projection
~~~~~~~~~~~~~~~~~~~~~
This script creates a light cone projection, a synthetic observation 
that stacks together projections from multiple datasets to extend over 
a given redshift interval.

.. yt_cookbook:: light_cone_projection.py

Light Cone with Halo Mask
~~~~~~~~~~~~~~~~~~~~~~~~~
This script combines the light cone generator with the halo profiler to 
make a light cone projection with all of the halos cut out of the image.

.. yt_cookbook:: light_cone_with_halo_mask.py 

Making Unique Light Cones
~~~~~~~~~~~~~~~~~~~~~~~~~
This script demonstrates how to make a series of light cone projections
that only have a maximum amount of volume in common.

.. yt_cookbook:: unique_light_cone_projections.py 

Making Light Rays
~~~~~~~~~~~~~~~~~
This script demonstrates how to make a synthetic quasar sight line and 
uses the halo profiler to record information about halos close to the 
line of sight.

.. yt_cookbook:: make_light_ray.py 

Creating and Fitting Absorption Spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This script demonstrates how to use light rays to create corresponding
absorption spectra and then fit the spectra to find absorbing
structures.

.. yt_cookbook:: fit_spectrum.py

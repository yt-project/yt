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

.. _cookbook-light_cone:

Light Cone Projection
~~~~~~~~~~~~~~~~~~~~~
This script creates a light cone projection, a synthetic observation 
that stacks together projections from multiple datasets to extend over 
a given redshift interval.

.. yt_cookbook:: light_cone_projection.py

.. _cookbook-light_ray:

Light Ray
~~~~~~~~~
This script demonstrates how to make a synthetic quasar sight line that 
extends over multiple datasets and can be used to generate a synthetic 
absorption spectrum.

.. yt_cookbook:: light_ray.py 

Creating and Fitting Absorption Spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This script demonstrates how to use light rays to create corresponding
absorption spectra and then fit the spectra to find absorbing
structures.

.. yt_cookbook:: fit_spectrum.py

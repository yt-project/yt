Cosmological Analysis
---------------------

These scripts demonstrate some basic and more advanced analysis that can be 
performed on cosmological simulation datasets.  Most of the following 
recipes are derived from functionality in yt's :ref:`analysis-modules`.

Plotting Halos
~~~~~~~~~~~~~~

This is a mechanism for plotting circles representing identified particle halos
on an image.
See :ref:`halo-analysis` and :ref:`annotate-halos` for more information.

.. yt_cookbook:: halo_plotting.py

.. _cookbook-rockstar-nested-grid:

Running Rockstar to Find Halos on Multi-Resolution-Particle Datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The version of Rockstar installed with yt does not have the capability
to work on datasets with particles of different masses.  Unfortunately,
many simulations possess particles of different masses, notably cosmological 
zoom datasets.  This recipe uses Rockstar in two different ways to generate a 
HaloCatalog from the highest resolution dark matter particles (the ones 
inside the zoom region).  It then overlays some of those halos on a projection
as a demonstration.  See :ref:`halo-analysis` and :ref:`annotate-halos` for
more information.

.. yt_cookbook:: rockstar_nest.py

.. _cookbook-halo_finding:

Halo Profiling and Custom Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script demonstrates the use of the halo catalog to create radial
profiles for each halo in a cosmological dataset.
See :ref:`halo_catalog` for more information.

.. yt_cookbook:: halo_profiler.py

.. _cookbook-light_cone:

Light Cone Projection
~~~~~~~~~~~~~~~~~~~~~

This script creates a light cone projection, a synthetic observation 
that stacks together projections from multiple datasets to extend over 
a given redshift interval.
See :ref:`light-cone-generator` for more information.

.. yt_cookbook:: light_cone_projection.py

.. _cookbook-light_ray:

Light Ray
~~~~~~~~~

This script demonstrates how to make a synthetic quasar sight line that 
extends over multiple datasets and can be used to generate a synthetic 
absorption spectrum.
See :ref:`light-ray-generator` and :ref:`absorption_spectrum` for more information.

.. yt_cookbook:: light_ray.py 

Creating and Fitting Absorption Spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script demonstrates how to use light rays to create corresponding
absorption spectra and then fit the spectra to find absorbing
structures.
See :ref:`light-ray-generator` and :ref:`absorption_spectrum` for more information.

.. yt_cookbook:: fit_spectrum.py

.. _examining-data:

Loading and Examining Data
==========================

Nominally, one should just be able to run ``yt.load()`` on a dataset and start
computing; however, there may be additional notes associated with different
data formats as described below.  Furthermore, we provide methods for loading
data from unsupported data formats in :ref:`loading-numpy-array`,
:ref:`generic-particle-data`, and :ref:`loading-spherical-data`.  Lastly, if
you want to examine the raw data for your particular dataset, visit
:ref:`low-level-data-inspection`.

.. toctree::
   :maxdepth: 2

   loading_data
   generic_array_data
   generic_particle_data
   spherical_data
   low_level_inspection

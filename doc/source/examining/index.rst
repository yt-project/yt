.. _examining-data:

Loading and Examining Data
==========================

Nominally, one should just be able to run ``yt.load()`` on a dataset and start
computing; however, there may be additional notes associated with different
data formats as described below.  Furthermore, we provide methods for loading
data from unsupported data formats in :doc:`Loading_Generic_Array_Data`,
:doc:`Loading_Generic_Particle_Data`, and :doc:`Loading_Spherical_Data`.
Lastly, if you want to examine the raw data for your particular dataset, visit
:ref:`low-level-data-inspection`.

.. toctree::
   :maxdepth: 2

   loading_data
   Loading_Generic_Array_Data
   Loading_Generic_Particle_Data
   Loading_Data_via_Functions
   Loading_Spherical_Data
   low_level_inspection

.. _analyzing:

General Data Analysis
=====================

This documentation describes much of the yt infrastructure for manipulating
one's data to extract the relevant information.  Fields, data objects, and
units are at the heart of how yt represents data.  Beyond this, we provide
a full description for how to filter your datasets based on specific criteria,
how to analyze chronological datasets from the same underlying simulation or
source (i.e. time series analysis), and how to run yt in parallel on
multiple processors to accomplish tasks faster.

.. toctree::
   :maxdepth: 2

   fields
   ../developing/creating_derived_fields
   objects
   units
   filtering
   generating_processed_data
   saving_data
   time_series_analysis
   particle_trajectories
   parallel_computation
   astropy_integrations

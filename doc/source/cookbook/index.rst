Cookbook
========

yt scripts can be a bit intimidating, and at times a bit obtuse.  But there's a
lot you can do, and this section of the manual will assist with figuring out
how to do some fairly common tasks -- which can lead to combining these, with
other Python code, into more complicated and advanced tasks.

There are some convenience functions to make writing standalone scripts a bit
easier.  Specifically, if you start your scripts with

.. code-block:: python

   from yt.mods import *
   pf = get_pf()

then your namespace is populated, and the last argument on the command line is
automatically turned into an :class:`EnzoStaticOutput`.

.. note::
   All of these scripts are located in the examples/ directory of the main
   distribution.  Except as otherwise noted, they are all executable if the last
   argument on the command line is the location of your parameter file:

.. code-block:: bash

   $ python2.5 examples/cookbook_hop_mass_sum.py /data/me/my_data0001

.. toctree::
   :maxdepth: 2

   making_simple_plots
   analyzing_data
   making_profiles
   plotting_profiles
   making_movies
   making_timeseries_plots
   extracting_data

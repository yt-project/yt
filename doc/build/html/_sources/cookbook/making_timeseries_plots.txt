Making Timeseries Plots
=======================

I don't think anyone will argue that timeseries data is important.  It is
possible to create it with yt; however, a few things should be noted.

Python is notorious for leaking memory via circular references -- if you
have a variable *var1* that references variable *var2*, when one is deleted, is
the other?  Efforts have been made to identify and remove all leaking
references, but if references are created in your code that generates the
timeseries data, you may find that more memory is used than is desired.

The process of generating time series data is fairly simple: you iterate over
the datasets, and generate a value for each.
(``cookbook_timeseries_max_dens.py``)

.. literalinclude:: ../../../examples/cookbook_timeseries_max_dens.py
   :language: python
   :linenos:

You could also imagine doing a time-evolution of an averaged quantity.  This
next example will (in a memory-conservative fashion) iterate over a bunch of
outputs and then return to you the averaged Temperature, weighted by Cell Mass.
(``cookbook_timeseries_avg_temp.py``)


.. literalinclude:: ../../../examples/cookbook_timeseries_avg_temp.py
   :language: python
   :linenos:

At this point, you can also write out the data to a file and plot it in another
plotting package:

.. code-block:: python
   :linenos:

   f = open("my_data.dat", "w")
   for i in range(len(times)):
       f.write("%0.5e %0.5e\n" % (times[i], avg_T[i]))
   f.close()

More complicated analysis can be done as well; :class:`PlotCollections`
created, slices made, etc.

The package :mod:`fido` inside of yt is also available, but it is not
well-documented.  It will store information about outputs, associated with a
given run, and you can iterate over those.  This is altogether a simpler
process, but it is not resilient against data movement.

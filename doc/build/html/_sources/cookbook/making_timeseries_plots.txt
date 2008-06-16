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
the datasets, and generate a value for each. ::

   >>> max_rho = []
   >>> max_pos = []
   >>> times = []
   >>> for i in range(30):
   ...     pf = lagos.EnzoStaticOutput("my_output%04i" % (i))
   ...     v, c = pf.h.find_max("Density")
   ...     max_rho.append(v)
   ...     max_pos.append(c)
   ...     times.append(pf["CosmologyCurrentTime"] * pf["years"])
   >>> pylab.loglog(times, max_rho)
   >>> pylab.clf()

You could also imagine doing a time-evolution of an averaged quantity::

   >>> avg_T = []
   >>> times = []
   >>> for i in range(30):
   ...     pf = lagos.EnzoStaticOutput("my_output%04i" % (i))
   ...     v, c = pf.h.find_max("Density")
   ...     sp = pf.h.sphere(c, 10.0/pf['kpc'])
   ...     avg_T.append(sp.quantities["WeightedAverageQuantity"]("Temperature",
   ...                                                           "CellMassMsun")
   ...     times.append(pf["CosmologyCurrentTime"] * pf["years"])
   >>> pylab.loglog(times, avg_T)

At this point, you can also write out the data to a file and plot it in another
plotting package: ::

   >>> f = open("my_data.dat", "w")
   >>> for i in range(len(times)):
   ...     f.write("%0.5e %0.5e\n" % (times[i], avg_T[i]))

More complicated analysis can be done as well; :class:`PlotCollections`
created, slices made, etc.

The package :mod:`fido` inside of yt is also available, but it is not
well-documented.  It will store information about outputs, associated with a
given run, and you can iterate over those.  This is altogether a simpler
process, but it is not resilient against data movement.

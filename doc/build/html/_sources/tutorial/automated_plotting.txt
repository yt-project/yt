.. index:: Plotting

Automated Plotting
==================

So now we know how to get and manipulate most of the datatypes, but we
still have to go through some processes to plot them manually.  However, yt has
facilities for plotting slices, projections and profiles in much easier ways.

.. index:: Raven, PlotCollections, Plotting;PlotCollections

Collections of Plots
--------------------

The plotting tool is called 'raven' and it's organized around the idea of
PlotCollection objects.  Each PlotCollection object is in some way thematically
associated -- maybe they are all the same data, sliced along different axes, for instance.
When you act on a PlotCollection, you act on all the plots that are part of that
collection.  Each PlotCollection is associated with a single parameter file.  ::

   >>> import yt.raven as raven
   >>> pc = raven.PlotCollection(a)

.. index:: Plotting;Adding plots to PlotCollections, Raven, PlotCollections;Adding Plots

Adding Plots
------------

Now we've got a plot collection -- a container for more plots that we can add
and then manipulate as a group.  So, let's try adding some plots to it!::

   >>> pc.add_slice("Density", 0)
   >>> pc.add_slice("Density", 1)
   >>> pc.add_slice("Density", 2)
   >>> pc.set_width(15, 'kpc')
   >>> pc.save("somename_15kpc")

There's a lot in that little block of text.  First I add a slice (note that while
we can feed it a pre-existing slice, it will also grab one automatically if we
don't!) along each axis (0,1,2) in the field "Density".  Each time we
add a plot to the PlotCollection, it is accessible as pc.plots[index] where 'index'
is the 0-based index in order of addition.  Each plot object has a 'data' property,
so if we wanted the EnzoSlice object, we can access it: ::

   >>> print pc.plots[1].data

After we've created the plots above, we set the width to 15 kiloparsecs.  Note
that raven understands any length unit that the hierarchy has, so you can set to
mpc, kpc, pc, au, km or cm.  Then we call save, and feed it a prefix -- it takes
care of adding some more information to the filename, so don't supply it a file
suffix. pc.save defaults to '.png', but a "format" keyword can be supplied::

   >>> pc.save("somename_15kpc",format='ps')

to get postscript output suitable for publication.

More Plots
----------

We have access to most of the interesting plots from the PlotCollection
interface, including cutting planes and projections.  (More on those later!)::

   >>> pc.add_projection("Temperature", 0, weight_field="CellMass")
   >>> pc.add_phase_sphere(1.0, 'pc', ["Density","Temperature","H2I_Fraction"])

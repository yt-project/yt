.. _manual-plotting:

Using the Manual Plotting Interface
===================================

Sometimes you need a lot of flexibility in creating plots. While the
:class:`~yt.visualization.plot_window.PlotWindow` provides an easy to
use object that can create nice looking, publication quality plots with a
minimum of effort, there are often times when its ease of use conflicts with
your need to change the font only on the x-axis, or whatever your
need/desire/annoying coauthor requires. To that end, ``yt`` provides a number of
ways of getting the raw data that goes into a plot to you in the form of a one
or two dimensional dataset that you can plot using any plotting method you like.
matplotlib or another python library are easiest, but these methods allow you to
take your data and plot it in gnuplot, or any unnamed commercial plotting
packages.

Note that the index object associated with your snapshot file contains a
list of plots you've made in ``ds.plots``.

.. _fixed-resolution-buffers:

Slice, Projections, and other Images: The Fixed Resolution Buffer
-----------------------------------------------------------------

For slices and projects, ``yt`` provides a manual plotting interface based on
the :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer` (hereafter
referred to as FRB) object. Despite its somewhat unwieldy name, at its heart, an
FRB is a very simple object: it's essentially a window into your data: you give
it a center and a width or a left and right edge, and an image resolution, and
the FRB returns a fully pixelized image. The simplest way to
generate an FRB is to use the ``.to_frb(width, resolution, center=None)`` method
of any data two-dimensional data object:

.. python-script::
   
   import pylab as P
   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

   c = pf.h.find_max('density')[1]
   proj = pf.proj('density', 0)

   width = (10, 'kpc') # we want a 1.5 mpc view
   res = [1000, 1000] # create an image with 1000x1000 pixels
   frb = proj.to_frb(width, res, center=c)

   P.imshow(np.array(frb['density']))
   P.savefig('my_perfect_figure.png')
   
The FRB is a very small object that can be deleted and recreated quickly (in
fact, this is how the reason GUI works when you pan and scan). Furthermore, you
can add new fields in the same "window", and each of them can be plotted with
their own zlimit. This is quite useful for creating a mosaic of the same region
in space with Density, Temperature, and x-velocity, for example. Each of these
quantities requires a substantially different set of limits.

A more complex example, showing a few ``yt`` helper functions that can make
setting up multiple axes with colorbars easier than it would be using only
matplotlib can be found in the cookbook.

.. _manual-line-plots:

Line Plots
----------

This is perhaps the simplest thing to do. ``yt`` provides a number of one dimensional objects, and these return a 1-D numpy array of their contents with direct dictionary access. As a simple example, take a :class:`~yt.data_objects.data_containers.AMROrthoRayBase` object, which can be created from a index by calling ``pf.ortho_ray(axis, center)``. 

.. python-script::

   from yt.mods import *
   import pylab as P
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   c = pf.h.find_max("density")[1]
   ax = 0 # take a line cut along the x axis
   ray = pf.ortho_ray(ax, (c[1], c[2])) # cutting through the y0,z0 such that we hit the max density

   P.subplot(211)
   P.semilogy(np.array(ray['x']), np.array(ray['density']))
   P.ylabel('density')
   P.subplot(212)
   P.semilogy(np.array(ray['x']), np.array(ray['temperature']))
   P.xlabel('x')
   P.ylabel('temperature')

   P.savefig("den_temp_xsweep.png")

Of course, you'll likely want to do something more sophisticated than using the
matplotlib defaults, but this gives the general idea.

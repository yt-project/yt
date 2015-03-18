.. _callbacks:

Plot Modifications: Overplotting Contours, Velocities, Particles, and More
==========================================================================

Adding callbacks to plots
-------------------------

Because the plots in yt are considered to be "volatile" -- existing
independent of the canvas on which they are plotted -- before they are saved,
you can have a set of "callbacks" run that modify them before saving to disk.
By adding a callback, you are telling the plot that whatever it does it itself,
your callback gets the last word.

Callbacks can be applied to plots created with
:class:`~yt.visualization.plot_window.SlicePlot`,
:class:`~yt.visualization.plot_window.ProjectionPlot`,
:class:`~yt.visualization.plot_window.OffAxisSlicePlot`, or
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot` by calling
one of the ``annotate_`` methods that hang off of the plot object.
The ``annotate_`` methods are dynamically generated based on the list
of available callbacks.  For example:

.. code-block:: python

   slc = SlicePlot(ds,0,'density')
   slc.annotate_title('This is a Density plot')

would add the :func:`~yt.visualization.plot_modifications.TitleCallback` to 
the plot object.  All of the callbacks listed below are available via 
similar ``annotate_`` functions.

To clear one or more annotations from an existing plot, see the 
:ref:`annotate_clear() function <annotate-clear>`.

Coordinate Systems in Callbacks
-------------------------------

Many of the callbacks (e.g. 
:class:`~yt.visualization.plot_modifications.TextLabelCallback`) are specified 
to occur at user-defined coordinate locations (like where to place a marker 
or text on the plot).  There are several different coordinate systems used 
to identify these locations.  These coordinate systems can be specified with 
the `coord_system` keyword in the relevant callback, which is by default 
set to `data`.  The valid coordinate systems are:

    `data` – the 3D dataset coordinates 
    `plot` – the 2D coordinates defined by the actual plot limits 
    `axis` – the MPL axis coordinates: (0,0) is lower left; (1,1) is upper right
    `figure` – the MPL figure coordinates: (0,0) is lower left, (1,1) is upper right

Here we will demonstrate these different coordinate systems for an projection 
of the x-plane (i.e. with axes in the y and z directions):

.. python-script::

import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
s = yt.SlicePlot(ds, 'x', 'density')
s.set_axes_unit('kpc')

# Plot marker and text in data coords
s.annotate_marker((0.2, 0.5, 0.9), coord_system='data')
s.annotate_text((0.2, 0.5, 0.9), 'data: (0.2, 0.5, 0.9)', coord_system='data')

# Plot marker and text in plot coords
s.annotate_marker((200, -300), coord_system='plot')
s.annotate_text((200, -300), 'plot: (200, -300)', coord_system='plot')

# Plot marker and text in axis coords
s.annotate_marker((0.1, 0.2), coord_system='axis')
s.annotate_text((0.1, 0.2), 'axis: (0.1, 0.2)', coord_system='axis')

# Plot marker and text in figure coords
# N.B. marker will not render outside of axis bounds
s.annotate_marker((0.1, 0.2), coord_system='figure', 
                  plot_args={'color':'black'})
s.annotate_text((0.1, 0.2), 'figure: (0.1, 0.2)', coord_system='figure', 
                text_args={'color':'black'})
s.save()

Available Callbacks
-------------------

The underlying functions are more thoroughly documented in :ref:`callback-api`.

.. include:: _cb_docstrings.inc

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

Available Callbacks
-------------------

The underlying functions are documented (largely identical to this) in
:ref:`callback-api`.

.. include:: _cb_docstrings.inc

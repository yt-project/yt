Primary Plot Interface
======================

Plot Collection
---------------

PlotCollection is the basic means by which most of your backend-plotting will
take place, and it contains a number of convenience functions for generating
images and manipulating existing plots.

.. currentmodule:: yt.raven

.. autoclass:: yt.raven.PlotCollection
   :members:

Plot Modification Callbacks
---------------------------

These are all meant to be instantiated and fed into
:meth:`yt.raven.RavenPlot.add_callback`.

.. autoclass:: yt.raven.PlotCallback

.. autoclass:: yt.raven.QuiverCallback

.. autoclass:: yt.raven.ParticleCallback

.. autoclass:: yt.raven.ContourCallback

.. autoclass:: yt.raven.GridBoundaryCallback

.. autoclass:: yt.raven.UnitBoundaryCallback

.. autoclass:: yt.raven.LinePlotCallback

.. autoclass:: yt.raven.CuttingQuiverCallback

3D Plots
--------

All of these plots require the usage of S2PLOT and its Python bindings.  Note,
also, that these require interaction with the user and do not support
disconnected use.  Furthermore, typically after a plot has been created the UI
remains onscreen for the remainder of the Python session.



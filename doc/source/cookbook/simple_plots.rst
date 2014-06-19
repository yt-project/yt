Making Simple Plots
-------------------

One of the easiest ways to interact with yt is by creating simple
visualizations of your data.  Below we show how to do this, as well as how to
extend these plots to be ready for publication.

Simple Slices
~~~~~~~~~~~~~

This script shows the simplest way to make a slice through a dataset.

.. yt_cookbook:: simple_slice.py

Simple Projections
~~~~~~~~~~~~~~~~~~

This is the simplest way to make a projection through a dataset.

.. yt_cookbook:: simple_projection.py

Simple Phase Plots
~~~~~~~~~~~~~~~~~~

This demonstrates how to make a phase plot.  Phase plots can be thought of as
two-dimensional histograms, where the value is either the weighted-average or
the total accumulation in a cell.

.. yt_cookbook:: simple_phase.py

Simple Probability Distribution Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often, one wants to examine the distribution of one variable as a function of
another.  This shows how to see the distribution of mass in a simulation, with
respect to the total mass in the simulation.

.. yt_cookbook:: simple_pdf.py

Simple 1D Histograms (Profiles)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a "profile," which is a 1D histogram.  This can be thought of as either
the total accumulation (when weight_field is set to ``None``) or the average 
(when a weight_field is supplied.)

.. yt_cookbook:: simple_profile.py

Simple Radial Profiles
~~~~~~~~~~~~~~~~~~~~~~

This shows how to make a profile of a quantity with respect to the radius, in
this case the radius in Mpc.

.. yt_cookbook:: simple_radial_profile.py

1D Profiles Over Time
~~~~~~~~~~~~~~~~~~~~~

This is a simple example of overplotting multiple 1D profiles from a number 
of datasets to show how they evolve over time.

.. yt_cookbook:: time_series_profiles.py

.. _cookbook-profile-variance:

Profiles with Variance Values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows how to plot the variance for a 1D profile.  In this example, we 
manually create a 1D profile object, which gives us access to the variance 
data.

.. yt_cookbook:: profile_with_variance.py

Making Plots of Multiple Fields Simultaneously
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By adding multiple fields to a single
:class:`~yt.visualization.plot_window.SlicePlot` or
:class:`~yt.visualization.plot_window.ProjectionPlot` some of the overhead of
creating the data object can be reduced, and better performance squeezed out.
This recipe shows how to add multiple fields to a single plot.

.. yt_cookbook:: simple_slice_with_multiple_fields.py 

Off-Axis Slicing
~~~~~~~~~~~~~~~~

A cutting plane allows you to slice at some angle that isn't aligned with the
axes.

.. yt_cookbook:: aligned_cutting_plane.py

.. _cookbook-simple-off-axis-projection:

Off-Axis Projection
~~~~~~~~~~~~~~~~~~~

Like cutting planes, off-axis projections can be created from any arbitrary 
viewing angle.

.. yt_cookbook:: simple_off_axis_projection.py

Simple Volume Rendering
~~~~~~~~~~~~~~~~~~~~~~~

Volume renderings are 3D projections rendering isocontours in any arbitrary
field (e.g. density, temperature, pressure, etc.)

.. yt_cookbook:: simple_volume_rendering.py

Showing and Hiding Axes Labels and Colorbars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example illustrates how to create a SlicePlot and then suppress the axes
labels and colorbars.  This is useful when you don't care about the physical
scales and just want to take a closer look at the raw plot data.

.. yt_cookbook:: show_hide_axes_colorbar.py

Accessing and Modifying Plots Directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While often the Plot Window, and its affiliated :ref:`callbacks` can
cover normal use cases, sometimes more direct access to the underlying
Matplotlib engine is necessary.  This recipe shows how to modify the plot
window :class:`matplotlib.axes.Axes` object directly.

.. yt_cookbook:: simple_slice_matplotlib_example.py 

.. _cookbook-simple_volume_rendering:

Image Background Colors
~~~~~~~~~~~~~~~~~~~~~~~

Here we see how to take an image and save it using different background colors. 

.. yt_cookbook:: image_background_colors.py

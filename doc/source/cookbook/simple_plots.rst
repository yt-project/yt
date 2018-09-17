Making Simple Plots
-------------------

One of the easiest ways to interact with yt is by creating simple
visualizations of your data.  Below we show how to do this, as well as how to
extend these plots to be ready for publication.

Simple Slices
~~~~~~~~~~~~~

This script shows the simplest way to make a slice through a dataset.  See
:ref:`slice-plots` for more information.

.. yt_cookbook:: simple_slice.py

Simple Projections (Non-Weighted)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the simplest way to make a projection through a dataset.  There are
several different :ref:`projection-types`, but non-weighted line integrals
and weighted line integrals are the two most common.  Here we create
density projections (non-weighted line integral).
See :ref:`projection-plots` for more information.

.. yt_cookbook:: simple_projection.py

Simple Projections (Weighted)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

And here we produce density-weighted temperature projections (weighted line
integral) for the same dataset as the non-weighted projections above.
See :ref:`projection-plots` for more information.

.. yt_cookbook:: simple_projection_weighted.py

Simple Phase Plots
~~~~~~~~~~~~~~~~~~

This demonstrates how to make a phase plot.  Phase plots can be thought of as
two-dimensional histograms, where the value is either the weighted-average or
the total accumulation in a cell.
See :ref:`how-to-make-2d-profiles` for more information.

.. yt_cookbook:: simple_phase.py

Simple 1D Line Plotting
~~~~~~~~~~~~~~~~~~~~~~~

This script shows how to make a ``LinePlot`` through a dataset.
See :ref:`manual-line-plots` for more information.

.. yt_cookbook:: simple_1d_line_plot.py

Simple Probability Distribution Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often, one wants to examine the distribution of one variable as a function of
another.  This shows how to see the distribution of mass in a simulation, with
respect to the total mass in the simulation.
See :ref:`how-to-make-2d-profiles` for more information.

.. yt_cookbook:: simple_pdf.py

Simple 1D Histograms (Profiles)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a "profile," which is a 1D histogram.  This can be thought of as either
the total accumulation (when weight_field is set to ``None``) or the average
(when a weight_field is supplied.)
See :ref:`how-to-make-1d-profiles` for more information.

.. yt_cookbook:: simple_profile.py

Simple Radial Profiles
~~~~~~~~~~~~~~~~~~~~~~

This shows how to make a profile of a quantity with respect to the radius.
See :ref:`how-to-make-1d-profiles` for more information.

.. yt_cookbook:: simple_radial_profile.py

1D Profiles Over Time
~~~~~~~~~~~~~~~~~~~~~

This is a simple example of overplotting multiple 1D profiles from a number
of datasets to show how they evolve over time.
See :ref:`how-to-make-1d-profiles` for more information.

.. yt_cookbook:: time_series_profiles.py

.. _cookbook-profile-stddev:

Profiles with Standard Deviation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows how to plot a 1D profile with error bars indicating the standard
deviation of the field values in each profile bin.  In this example, we manually
create a 1D profile object, which gives us access to the standard deviation
data.  See :ref:`how-to-make-1d-profiles` for more information.

.. yt_cookbook:: profile_with_standard_deviation.py

Making Plots of Multiple Fields Simultaneously
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By adding multiple fields to a single
:class:`~yt.visualization.plot_window.SlicePlot` or
:class:`~yt.visualization.plot_window.ProjectionPlot` some of the overhead of
creating the data object can be reduced, and better performance squeezed out.
This recipe shows how to add multiple fields to a single plot.
See :ref:`slice-plots` and :ref:`projection-plots` for more information.

.. yt_cookbook:: simple_slice_with_multiple_fields.py

Off-Axis Slicing
~~~~~~~~~~~~~~~~

One can create slices from any arbitrary angle, not just those aligned with
the x,y,z axes.
See :ref:`off-axis-slices` for more information.

.. yt_cookbook:: simple_off_axis_slice.py

.. _cookbook-simple-off-axis-projection:

Off-Axis Projection
~~~~~~~~~~~~~~~~~~~

Like off-axis slices, off-axis projections can be created from any arbitrary
viewing angle.
See :ref:`off-axis-projections` for more information.

.. yt_cookbook:: simple_off_axis_projection.py

.. _cookbook-simple-particle-plot:

Simple Particle Plot
~~~~~~~~~~~~~~~~~~~~

You can also use yt to make particle-only plots. This script shows how to
plot all the particle x and y positions in a dataset, using the particle mass
to set the color scale.
See :ref:`particle-plots` for more information.

.. yt_cookbook:: particle_xy_plot.py

.. _cookbook-non-spatial-particle-plot:

Non-spatial Particle Plots
~~~~~~~~~~~~~~~~~~~~~~~~~~

You are not limited to plotting spatial fields on the x and y axes. This
example shows how to plot the particle x-coordinates versus their z-velocities,
again using the particle mass to set the colorbar.
See :ref:`particle-plots` for more information.

.. yt_cookbook:: particle_xvz_plot.py

.. _cookbook-single-color-particle-plot:

Single-color Particle Plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't want to display a third field on the color bar axis, simply pass
in a color string instead of a particle field.
See :ref:`particle-plots` for more information.

.. yt_cookbook:: particle_one_color_plot.py

.. _cookbook-simple_volume_rendering:

Simple Volume Rendering
~~~~~~~~~~~~~~~~~~~~~~~

Volume renderings are 3D projections rendering isocontours in any arbitrary
field (e.g. density, temperature, pressure, etc.)
See :ref:`volume_rendering` for more information.

.. yt_cookbook:: simple_volume_rendering.py

.. _show-hide-axes-colorbar:

Showing and Hiding Axis Labels and Colorbars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example illustrates how to create a SlicePlot and then suppress the axes
labels and colorbars.  This is useful when you don't care about the physical
scales and just want to take a closer look at the raw plot data.  See
:ref:`hiding-colorbar-and-axes` for more information.

.. yt_cookbook:: show_hide_axes_colorbar.py


.. _cookbook_label_formats:

Setting Field Label Formats
---------------------------

This example illustrates how to change the label format for
ion species from the default roman numeral style.

.. yt_cookbook:: changing_label_formats.py


.. _matplotlib-primitives:

Accessing and Modifying Plots Directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While often the Plot Window, and its affiliated :ref:`callbacks` can
cover normal use cases, sometimes more direct access to the underlying
Matplotlib engine is necessary.  This recipe shows how to modify the plot
window :class:`matplotlib.axes.Axes` object directly.
See :ref:`matplotlib-customization` for more information.

.. yt_cookbook:: simple_slice_matplotlib_example.py

Changing the Colormap used in a Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yt has sensible defaults for colormaps, but there are over a hundred available
for customizing your plots.  Here we generate a projection and then change
its colormap.  See :ref:`colormaps` for a list and for images of all the
available colormaps.

.. yt_cookbook:: colormaps.py

Image Background Colors
~~~~~~~~~~~~~~~~~~~~~~~

Here we see how to take an image and save it using different background colors.

In this case we use the :ref:`cookbook-simple_volume_rendering`
recipe to generate the image, but it works for any NxNx4 image array
(3 colors and 1 opacity channel).  See :ref:`volume_rendering` for more
information.

.. yt_cookbook:: image_background_colors.py

.. _annotations-recipe:

Annotating Plots to Include Lines, Text, Shapes, etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can be useful to add annotations to plots to show off certain features
and make it easier for your audience to understand the plot's purpose.  There
are a variety of available :ref:`plot modifications <callbacks>` one can use
to add annotations to their plots.  Below includes just a handful, but please
look at the other :ref:`plot modifications <callbacks>` to get a full
description of what you can do to highlight your figures.

.. yt_cookbook:: annotations.py

Annotating Plots with a Timestamp and Physical Scale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When creating movies of multiple outputs from the same simulation (see :ref:`time-series-analysis`), it can be helpful to include a timestamp and the physical scale of each individual output.  This is simply achieved using the :ref:`annotate_timestamp() <annotate-timestamp>` and :ref:`annotate_scale() <annotate-scale>` callbacks on your plots.  For more information about similar plot modifications using other callbacks, see the section on :ref:`Plot Modifications <callbacks>`.

.. yt_cookbook:: annotate_timestamp_and_scale.py

A Few Complex Plots
-------------------

The built-in plotting functionality covers the very simple use cases that are
most common.  These scripts will demonstrate how to construct more complex
plots or publication-quality plots.  In many cases these show how to make
multi-panel plots.

Multi-Width Image
~~~~~~~~~~~~~~~~~

This is a simple recipe to show how to open a dataset and then plot slices
through it at varying widths.
See :ref:`slice-plots` for more information.

.. yt_cookbook:: multi_width_image.py

.. _image-resolution-primer:

Varying the resolution of an image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This illustrates the various parameters that control the resolution
of an image, including the (deprecated) refinement level, the size of
the :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`,
and the number of pixels in the output image.

In brief, there are three parameters that control the final resolution,
with a fourth entering for particle data that is deposited onto a mesh
(i.e. pre-4.0).  Those are:

1. ``buff_size``, which can be altered with
:meth:`~yt.visualization.plot_window.PlotWindow.set_buff_size`, which
is inherited by
:class:`~yt.visualization.plot_window.AxisAlignedSlicePlot`,
:class:`~yt.visualization.plot_window.OffAxisSlicePlot`,
:class:`~yt.visualization.plot_window.ProjectionPlot`, and
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`.  This
controls the number of resolution elements in the
:class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`,
which can be thought of as the number of individually colored
squares (on a side) in a 2D image. ``buff_size`` can be set
after creating the image with
:meth:`~yt.visualization.plot_window.PlotWindow.set_buff_size`,
or during image creation with the ``buff_size`` argument to any
of the four preceding classes.

2. ``figure_size``, which can be altered with either
:meth:`~yt.visualization.plot_container.PlotContainer.set_figure_size`
or with :meth:`~yt.visualization.plot_container.PlotWindow.set_window_size`
(the latter simply calls
:meth:`~yt.visualization.plot_container.PlotContainer.set_figure_size`),
or can be set during image creation with the ``window_size`` argument.
This sets the size of the final image (including the visualization and,
if applicable, the axes and colorbar as well) in inches.

3. ``dpi``, i.e. the dots-per-inch in your final file, which can also
be thought of as the actual resolution of your image.  This can
only be set on save via the ``mpl_kwargs`` parameter to
:meth:`~yt.visualization.plot_container.PlotContainer.save`.  The
``dpi`` and ``figure_size`` together set the true resolution of your
image (final image will be ``dpi`` :math:`*` ``figure_size`` pixels on a
side), so if these are set too low, then your ``buff_size`` will not
matter.  On the other hand, increasing these without increasing
``buff_size`` accordingly will simply blow up your resolution
elements to fill several real pixels.

4. (only for meshed particle data) ``n_ref``, the maximum nubmer of
particles in a cell in the oct-tree allowed before it is refined
(removed in yt-4.0 as particle data is no longer deposited onto
an oct-tree).  For particle data, ``n_ref`` effectively sets the
underlying resolution of your simulation.  Regardless, for either
grid data or deposited particle data, your image will never be
higher resolution than your simulation data.  In other words,
if you are visualizing a region 50 kpc across that includes
data that reaches a resolution of 100 pc, then there's no reason
to set a ``buff_size`` (or a ``dpi`` :math:`*` ``figure_size``) above
50 kpc/ 100 pc = 500.

The below script demonstrates how each of these can be varied.

.. yt_cookbook:: image_resolution.py


Multipanel with Axes Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This illustrates how to use a SlicePlot to control a multipanel plot.  This
plot uses axes labels to illustrate the length scales in the plot.
See :ref:`slice-plots` and the
`Matplotlib AxesGrid Object <https://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html>`_
for more information.

.. yt_cookbook:: multiplot_2x2.py

The above example gives you full control over the plots, but for most
purposes, the ``export_to_mpl_figure`` method is a simpler option,
allowing us to make a similar plot as:

.. yt_cookbook:: multiplot_export_to_mpl.py

Multipanel with PhasePlot
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This illustrates how to use PhasePlot in a multipanel plot.
See :ref:`how-to-make-2d-profiles` and the
`Matplotlib AxesGrid Object <https://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html>`_
for more information.

.. yt_cookbook:: multiplot_phaseplot.py

Time Series Multipanel
~~~~~~~~~~~~~~~~~~~~~~

This illustrates how to create a multipanel plot of a time series dataset.
See :ref:`projection-plots`, :ref:`time-series-analysis`, and the
`Matplotlib AxesGrid Object <https://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html>`_
for more information.

.. yt_cookbook:: multiplot_2x2_time_series.py

Multiple Slice Multipanel
~~~~~~~~~~~~~~~~~~~~~~~~~

This illustrates how to create a multipanel plot of slices along the coordinate
axes.  To focus on what's happening in the x-y plane, we make an additional
Temperature slice for the bottom-right subpanel.
See :ref:`slice-plots` and the
`Matplotlib AxesGrid Object <https://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html>`_
for more information.

.. yt_cookbook:: multiplot_2x2_coordaxes_slice.py

Multi-Plot Slice and Projections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows how to combine multiple slices and projections into a single image,
with detailed control over colorbars, titles and color limits.
See :ref:`slice-plots` and :ref:`projection-plots` for more information.

.. yt_cookbook:: multi_plot_slice_and_proj.py

.. _advanced-multi-panel:

Advanced Multi-Plot Multi-Panel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This produces a series of slices of multiple fields with different color maps
and zlimits, and makes use of the FixedResolutionBuffer. While this is more
complex than the equivalent plot collection-based solution, it allows for a
*lot* more flexibility. Every part of the script uses matplotlib commands,
allowing its full power to be exercised.
See :ref:`slice-plots` and :ref:`projection-plots` for more information.

.. yt_cookbook:: multi_plot_3x2_FRB.py

Time Series Movie
~~~~~~~~~~~~~~~~~

This shows how to use matplotlib's animation framework with yt plots.

.. yt_cookbook:: matplotlib-animation.py

.. _cookbook-offaxis_projection:

Off-Axis Projection (an alternate method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to take an image-plane line integral along an
arbitrary axis in a simulation.  This uses alternate machinery than the
standard :ref:`PlotWindow interface <off-axis-projections>` to create an
off-axis projection as demonstrated in this
:ref:`recipe <cookbook-simple-off-axis-projection>`.

.. yt_cookbook:: offaxis_projection.py

Off-Axis Projection with a Colorbar (an alternate method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe shows how to generate a colorbar with a projection of a dataset
from an arbitrary projection angle (so you are not confined to the x, y, and z
axes).

This uses alternate machinery than the standard
:ref:`PlotWindow interface <off-axis-projections>` to create an off-axis
projection as demonstrated in this
:ref:`recipe <cookbook-simple-off-axis-projection>`.

.. yt_cookbook:: offaxis_projection_colorbar.py

.. _thin-slice-projections:

Thin-Slice Projections
~~~~~~~~~~~~~~~~~~~~~~

This recipe is an example of how to project through only a given data object,
in this case a thin region, and then display the result.
See :ref:`projection-plots` and :ref:`available-objects` for more information.

.. yt_cookbook:: thin_slice_projection.py

Plotting Particles Over Fluids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to overplot particles on top of a fluid image.
See :ref:`annotate-particles` for more information.

.. yt_cookbook:: overplot_particles.py

Plotting Grid Edges Over Fluids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to overplot grid boxes on top of a fluid image.
Each level is represented with a different color from white (low refinement) to
black (high refinement).  One can change the colormap used for the grids colors
by using the cmap keyword (or set it to None to get all grid edges as black).
See :ref:`annotate-grids` for more information.

.. yt_cookbook:: overplot_grids.py

Overplotting Velocity Vectors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to plot velocity vectors on top of a slice.
See :ref:`annotate-velocity` for more information.

.. yt_cookbook:: velocity_vectors_on_slice.py

Overplotting Contours
~~~~~~~~~~~~~~~~~~~~~

This is a simple recipe to show how to open a dataset, plot a slice through it,
and add contours of another quantity on top.
See :ref:`annotate-contours` for more information.

.. yt_cookbook:: contours_on_slice.py

Simple Contours in a Slice
~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it is useful to plot just a few contours of a quantity in a
dataset.  This shows how one does this by first making a slice, adding
contours, and then hiding the colormap plot of the slice to leave the
plot containing only the contours that one has added.
See :ref:`annotate-contours` for more information.

.. yt_cookbook:: simple_contour_in_slice.py

Styling Radial Profile Plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates a method of calculating radial profiles for several
quantities, styling them and saving out the resultant plot.
See :ref:`how-to-make-1d-profiles` for more information.

.. yt_cookbook:: radial_profile_styles.py

Customized Profile Plot
~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to create a fully customized 1D profile object
using the :func:`~yt.data_objects.profiles.create_profile` function and then
create a :class:`~yt.visualization.profile_plotter.ProfilePlot` using the
customized profile.  This illustrates how a ``ProfilePlot`` created this way
inherits the properties of the profile it is constructed from.
See :ref:`how-to-make-1d-profiles` for more information.

.. yt_cookbook:: customized_profile_plot.py

Customized Phase Plot
~~~~~~~~~~~~~~~~~~~~~

Similar to the recipe above, this demonstrates how to create a fully customized
2D profile object using the :func:`~yt.data_objects.profiles.create_profile`
function and then create a :class:`~yt.visualization.profile_plotter.PhasePlot`
using the customized profile object.  This illustrates how a ``PhasePlot``
created this way inherits the properties of the profile object from which it
is constructed. See :ref:`how-to-make-2d-profiles` for more information.

.. yt_cookbook:: customized_phase_plot.py

.. _cookbook-camera_movement:

Moving a Volume Rendering Camera
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this recipe, we move a camera through a domain and take multiple volume
rendering snapshots. This recipe uses an unstructured mesh dataset (see
:ref:`unstructured_mesh_rendering`), which makes it easier to visualize what
the Camera is doing, but you can manipulate the Camera for other dataset types
in exactly the same manner.

See :ref:`camera_movement` for more information.

.. yt_cookbook:: camera_movement.py

Volume Rendering with Custom Camera
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this recipe we modify the :ref:`cookbook-simple_volume_rendering` recipe to
use customized camera properties. See :ref:`volume_rendering` for more
information.

.. yt_cookbook:: custom_camera_volume_rendering.py

.. _cookbook-custom-transfer-function:

Volume Rendering with a Custom Transfer Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this recipe we modify the :ref:`cookbook-simple_volume_rendering` recipe to
use customized camera properties. See :ref:`volume_rendering` for more
information.

.. yt_cookbook:: custom_transfer_function_volume_rendering.py

.. _cookbook-sigma_clip:

Volume Rendering with Sigma Clipping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this recipe we output several images with different values of sigma_clip
set in order to change the contrast of the resulting image.  See
:ref:`sigma_clip` for more information.

.. yt_cookbook:: sigma_clip.py

Zooming into an Image
~~~~~~~~~~~~~~~~~~~~~

This is a recipe that takes a slice through the most dense point, then creates
a bunch of frames as it zooms in.  It's important to note that this particular
recipe is provided to show how to be more flexible and add annotations and the
like -- the base system, of a zoomin, is provided by the "yt zoomin" command on
the command line.
See :ref:`slice-plots` and :ref:`callbacks` for more information.

.. yt_cookbook:: zoomin_frames.py

.. _cookbook-various_lens:

Various Lens Types for Volume Rendering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example illustrates the usage and feature of different lenses for volume rendering.

.. yt_cookbook:: various_lens.py

.. _cookbook-opaque_rendering:

Opaque Volume Rendering
~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to make semi-opaque volume renderings, but also
how to step through and try different things to identify the type of volume
rendering you want.
See :ref:`opaque_rendering` for more information.

.. yt_cookbook:: opaque_rendering.py

Volume Rendering Multiple Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can render multiple fields by adding new ``VolumeSource`` objects to the
scene for each field you want to render.

.. yt_cookbook:: render_two_fields.py

.. _cookbook-amrkdtree_downsampling:

Downsampling Data for Volume Rendering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to downsample data in a simulation to speed up
volume rendering.
See :ref:`volume_rendering` for more information.

.. yt_cookbook:: amrkdtree_downsampling.py

.. _cookbook-volume_rendering_annotations:

Volume Rendering with Bounding Box and Overlaid Grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to overplot a bounding box on a volume rendering
as well as overplotting grids representing the level of refinement achieved
in different regions of the code.
See :ref:`volume_rendering_annotations` for more information.

.. yt_cookbook:: rendering_with_box_and_grids.py

Volume Rendering with Annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to write the simulation time, show an
axis triad indicating the direction of the coordinate system, and show
the transfer function on a volume rendering.  Please note that this
recipe relies on the old volume rendering interface.  While one can
continue to use this interface, it may be incompatible with some of the
new developments and the infrastructure described in :ref:`volume_rendering`.

.. yt_cookbook:: vol-annotated.py

.. _cookbook-vol-points:

Volume Rendering with Points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to make a volume rendering composited with point
sources. This could represent star or dark matter particles, for example.

.. yt_cookbook:: vol-points.py

.. _cookbook-vol-lines:

Volume Rendering with Lines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to make a volume rendering composited with line
sources.

.. yt_cookbook:: vol-lines.py

Plotting Streamlines
~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to display streamlines in a simulation.  (Note:
streamlines can also be queried for values!)
See :ref:`streamlines` for more information.

.. yt_cookbook:: streamlines.py

Plotting Isocontours
~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to extract an isocontour and then plot it in
matplotlib, coloring the surface by a second quantity.
See :ref:`surfaces` for more information.

.. yt_cookbook:: surface_plot.py

Plotting Isocontours and Streamlines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe plots both isocontours and streamlines simultaneously.  Note that
this will not include any blending, so streamlines that are occluded by the
surface will still be visible.
See :ref:`streamlines` and :ref:`surfaces` for more information.

.. yt_cookbook:: streamlines_isocontour.py

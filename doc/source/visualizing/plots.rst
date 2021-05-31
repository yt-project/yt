
.. _how-to-make-plots:

How to Make Plots
=================

In this section we explain how to use yt to create visualizations
of simulation data, derived fields, and the data produced by yt
analysis objects.  For details about the data extraction and
algorithms used to produce the image and analysis data, please see the
yt `method paper
<https://ui.adsabs.harvard.edu/abs/2011ApJS..192....9T>`_.  There are also
many example scripts in :ref:`cookbook`.

The :class:`~yt.visualization.plot_window.PlotWindow` interface is useful for
taking a quick look at simulation outputs.  Simple mechanisms exist for making
plots of slices, projections, 1D spatial line plots, 1D profiles, and 2D
profiles (phase plots), all of which are described below.

.. _viewing-plots:

Viewing Plots
-------------

YT uses an environment neutral plotting mechanism that detects the appropriate
matplotlib configuration for a given environment, however it defaults to a basic
renderer. To utilize interactive plots in matplotlib supported
environments (Qt, GTK, WX, etc.) simply call the ``toggle_interactivity()`` function. Below is an
example in a jupyter notebook environment, but the same command should work
in other environments as well:

.. code-block:: IPython

   %matplotlib notebook
   import yt
   yt.toggle_interactivity()

.. _simple-inspection:

Slices & Projections
--------------------

If you need to take a quick look at a single simulation output, yt
provides the :class:`~yt.visualization.plot_window.PlotWindow` interface for
generating annotated 2D visualizations of simulation data.  You can create a
:class:`~yt.visualization.plot_window.PlotWindow` plot by
supplying a dataset, a list of fields to plot, and a plot center to
create a :class:`~yt.visualization.plot_window.AxisAlignedSlicePlot`,
:class:`~yt.visualization.plot_window.OffAxisSlicePlot`,
:class:`~yt.visualization.plot_window.ProjectionPlot`, or
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`.

Plot objects use yt data objects to extract the maximum resolution
data available to render a 2D image of a field. Whenever a
two-dimensional image is created, the plotting object first obtains
the necessary data at the *highest resolution*.  Every time an image
is requested of it -- for instance, when the width or field is changed
-- this high-resolution data is then pixelized and placed in a buffer
of fixed size. This is accomplished behind the scenes using
:class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`.

The :class:`~yt.visualization.plot_window.PlotWindow` class exposes the
underlying matplotlib
`figure <https://matplotlib.org/stable/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure>`_
and `axes <https://matplotlib.org/stable/api/axes_api.html#matplotlib.axes.Axes>`_
objects, making it easy to customize your plots and
add new annotations.  See :ref:`matplotlib-customization` for more information.

.. _slice-plots:

Slice Plots
~~~~~~~~~~~

The quickest way to plot a slice of a field through your data is via
:class:`~yt.visualization.plot_window.SlicePlot`.  These plots are generally
quicker than projections because they only need to read and process a slice
through the dataset.

The following script plots a slice through the density field along the z-axis
centered on the center of the simulation box in a simulation dataset we've
opened and stored in ``ds``:

.. code-block:: python

    slc = yt.SlicePlot(ds, "z", "density")
    slc.save()

These two commands will create a slice object and store it in a variable we've
called ``slc``.  Since this plot is aligned with the simulation coordinate
system, ``slc`` is an instance of
:class:`~yt.visualization.plot_window.AxisAlignedSlicePlot`. We then call the
``save()`` function, which automatically saves the plot in png image format with
an automatically generated filename.  If you don't want the slice object to
stick around, you can accomplish the same thing in one line:

.. code-block:: python

    yt.SlicePlot(ds, "z", "density").save()

It's nice to keep the slice object around if you want to modify the plot.  By
default, the plot width will be set to the size of the simulation box.  To zoom
in by a factor of ten, you can call the zoom function attached to the slice
object:

.. code-block:: python

    slc = yt.SlicePlot(ds, "z", "density")
    slc.zoom(10)
    slc.save("zoom")

This will save a new plot to disk with a different filename - prepended with
'zoom' instead of the name of the dataset. If you want to set the width
manually, you can do that as well. For example, the following sequence of
commands will create a slice, set the width of the plot to 10 kiloparsecs, and
save it to disk.

.. code-block:: python

    from yt.units import kpc

    slc = yt.SlicePlot(ds, "z", "density")
    slc.set_width(10 * kpc)
    slc.save("10kpc")

The plot width can be specified independently along the x and y direction by
passing a tuple of widths.  An individual width can also be represented using a
``(value, unit)`` tuple.  The following sequence of commands all equivalently
set the width of the plot to 200 kiloparsecs in the ``x`` and ``y`` direction.

.. code-block:: python

    from yt.units import kpc

    slc.set_width(200 * kpc)
    slc.set_width((200, "kpc"))
    slc.set_width((200 * kpc, 200 * kpc))

The ``SlicePlot`` also optionally accepts the coordinate to center the plot on
and the width of the plot:

.. code-block:: python

    yt.SlicePlot(ds, "z", "density", center=[0.2, 0.3, 0.8], width=(10, "kpc")).save()

Note that, by default,
:class:`~yt.visualization.plot_window.SlicePlot` shifts the
coordinates on the axes such that the origin is at the center of the
slice.  To instead use the coordinates as defined in the dataset, use
the optional argument: ``origin="native"``

If supplied without units, the center is assumed by in code units.  There are also
the following alternative options for the ``center`` keyword:

* ``"center"``, ``"c"``: the domain center
* ``"max"``, ``"m"``: the position of the maximum density
* ``("min", field)``: the position of the minimum of ``field``
* ``("max", field)``: the position of the maximum of ``field``

where for the last two objects any spatial field, such as ``"density"``,
``"velocity_z"``,
etc., may be used, e.g. ``center=("min","temperature")``.

The effective resolution of the plot (i.e. the number of resolution elements
in the image itself) can be controlled with the ``buff_size`` argument:

.. code-block:: python

    yt.SlicePlot(ds, "z", "density", buff_size=(1000, 1000))


Here is an example that combines all of the options we just discussed.

.. python-script::

   import yt
   from yt.units import kpc

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(
       ds,
       "z",
       "density",
       center=[0.5, 0.5, 0.5],
       width=(20, "kpc"),
       buff_size=(1000, 1000),
   )
   slc.save()

The above example will display an annotated plot of a slice of the
Density field in a 20 kpc square window centered on the coordinate
(0.5, 0.5, 0.5) in the x-y plane.  The axis to slice along is keyed to the
letter 'z', corresponding to the z-axis.  Finally, the image is saved to
a png file.

Conceptually, you can think of the plot object as an adjustable window
into the data. For example:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "pressure", center="c")
   slc.save()
   slc.zoom(30)
   slc.save("zoom")

will save a plot of the pressure field in a slice along the z
axis across the entire simulation domain followed by another plot that
is zoomed in by a factor of 30 with respect to the original
image. Both plots will be centered on the center of the simulation box.
With these sorts of manipulations, one can easily pan and zoom onto an
interesting region in the simulation and adjust the boundaries of the
region to visualize on the fly.

If you want to slice through a subset of the full dataset volume,
you can use the ``data_source`` keyword with a :ref:`data object <data-objects>`
or a :ref:`cut region <cut-regions>`.

See :class:`~yt.visualization.plot_window.AxisAlignedSlicePlot` for the
full class description.

.. _plot-2d:

Plots of 2D Datasets
~~~~~~~~~~~~~~~~~~~~

If you have a two-dimensional cartesian, cylindrical, or polar dataset,
:func:`~yt.visualization.plot_window.plot_2d` is a way to make a plot
within the dataset's plane without having to specify the axis, which
in this case is redundant. Otherwise, ``plot_2d`` accepts the same
arguments as ``SlicePlot``. The one other difference is that the
``center`` keyword argument can be a two-dimensional coordinate instead
of a three-dimensional one:

.. python-script::

    import yt

    ds = yt.load("WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030")
    p = yt.plot_2d(ds, "density", center=[1.0, 0.4])
    p.set_log("density", False)
    p.save()

See :func:`~yt.visualization.plot_window.plot_2d` for the full description
of the function and its keywords.

.. _off-axis-slices:

Off Axis Slices
~~~~~~~~~~~~~~~

Off axis slice plots can be generated in much the same way as
grid-aligned slices.  Off axis slices use
:class:`~yt.data_objects.selection_data_containers.YTCuttingPlane` to slice
through simulation domains at an arbitrary oblique angle.  A
:class:`~yt.visualization.plot_window.OffAxisSlicePlot` can be
instantiated by specifying a dataset, the normal to the cutting
plane, and the name of the fields to plot.  Just like an
:class:`~yt.visualization.plot_window.AxisAlignedSlicePlot`, an
:class:`~yt.visualization.plot_window.OffAxisSlicePlot` can be created via the
:class:`~yt.visualization.plot_window.SlicePlot` class. For example:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1, 1, 0]  # vector normal to cutting plane
   north_vector = [-1, 1, 0]
   cut = yt.SlicePlot(ds, L, "density", width=(25, "kpc"), north_vector=north_vector)
   cut.save()

In this case, a normal vector for the cutting plane is supplied in the second
argument. Optionally, a ``north_vector`` can be specified to fix the orientation
of the image plane.

.. _projection-plots:

Projection Plots
~~~~~~~~~~~~~~~~

Using a fast adaptive projection, yt is able to quickly project
simulation data along the coordinate axes.

Projection plots are created by instantiating a
:class:`~yt.visualization.plot_window.ProjectionPlot` object.  For
example:

.. python-script::

   import yt
   from yt.units import kpc

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   prj = yt.ProjectionPlot(
       ds,
       2,
       "temperature",
       width=25 * kpc,
       weight_field="density",
       buff_size=(1000, 1000),
   )
   prj.save()

will create a density-weighted projection of the temperature field along
the x axis with 1000 resolution elements per side, plot it, and then save
the plot to a png image file.

Like :ref:`slice-plots`, annotations and modifications can be applied
after creating the ``ProjectionPlot`` object.  Annotations are
described in :ref:`callbacks`.  See
:class:`~yt.visualization.plot_window.ProjectionPlot` for the full
class description.

If you want to project through a subset of the full dataset volume,
you can use the ``data_source`` keyword with a :ref:`data object <data-objects>`.
The :ref:`thin-slice-projections` recipes demonstrates this functionality.

.. _projection-types:

Types of Projections
""""""""""""""""""""

There are several different methods of projections that can be made either
when creating a projection with ds.proj() or when making a ProjectionPlot.
In either construction method, set the ``method`` keyword to be one of the
following:

``integrate`` (unweighted)
    This is the default projection method. It simply integrates the
    requested field  :math:`f(x)` along a line of sight  :math:`\hat{n}` ,
    given by the axis parameter (e.g. :math:`\hat{i},\hat{j},` or
    :math:`\hat{k}`).  The units of the projected field
    :math:`g(X)` will be the units of the unprojected field  :math:`f(x)`
    multiplied by the appropriate length unit, e.g., density in
    :math:`\mathrm{g\ cm^{-3}}` will be projected to  :math:`\mathrm{g\ cm^{-2}}`.

.. math::

    g(X) = {\int\ {f(x)\hat{n}\cdot{dx}}}

``integrate`` (weighted)
    When using the ``integrate``  method, a ``weight_field`` argument may also
    be specified, which will produce a weighted projection.  :math:`w(x)`
    is the field used as a weight. One common example would
    be to weight the "temperature" field by the "density" field. In this case,
    the units of the projected field are the same as the unprojected field.

.. math::

    g(X) = \frac{\int\ {f(x)w(x)\hat{n}\cdot{dx}}}{\int\ {w(x)\hat{n}\cdot{dx}}}

``mip``
    This method picks out the maximum value of a field along the line of
    sight given by the axis parameter.

``sum``
    This method is the same as ``integrate``, except that it does not
    multiply by a path length when performing the integration, and is just a
    straight summation of the field along the given axis. The units of the
    projected field will be the same as those of the unprojected field. This
    method is typically only useful for datasets such as 3D FITS cubes where
    the third axis of the dataset is something like velocity or frequency, and
    should _only_ be used with fixed-resolution grid-based datasets.

.. _off-axis-projections:

Off Axis Projection Plots
~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, off axis projections are created using :ref:`camera`
by applying the
:class:`~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction`.
In this use case, the volume renderer casts a set of plane parallel rays, one
for each pixel in the image.  The data values along each ray are summed,
creating the final image buffer.

.. _off-axis-projection-function:

To avoid manually creating a camera and setting the transfer
function, yt provides the
:func:`~yt.visualization.volume_rendering.off_axis_projection.off_axis_projection`
function, which wraps the camera interface to create an off axis
projection image buffer.  These images can be saved to disk or
used in custom plots.  This snippet creates an off axis
projection through a simulation.

.. python-script::

   import numpy as np

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1, 1, 0]  # vector normal to cutting plane
   north_vector = [-1, 1, 0]
   W = [0.02, 0.02, 0.02]
   c = [0.5, 0.5, 0.5]
   N = 512
   image = yt.off_axis_projection(ds, c, L, W, N, "density")
   yt.write_image(np.log10(image), "%s_offaxis_projection.png" % ds)

Here, ``W`` is the width of the projection in the x, y, *and* z
directions.

One can also generate annotated off axis projections using
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`. These
plots can be created in much the same way as an
``OffAxisSlicePlot``, requiring only an open dataset, a direction
to project along, and a field to project.  For example:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1, 1, 0]  # vector normal to cutting plane
   north_vector = [-1, 1, 0]
   prj = yt.OffAxisProjectionPlot(
       ds, L, "density", width=(25, "kpc"), north_vector=north_vector
   )
   prj.save()

OffAxisProjectionPlots can also be created with a number of
keyword arguments, as described in
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`

.. _unstructured-mesh-slices:

Unstructured Mesh Slices
------------------------

Unstructured Mesh datasets can be sliced using the same syntax as above.
Here is an example script using a publicly available MOOSE dataset:

.. python-script::

   import yt

   ds = yt.load("MOOSE_sample_data/out.e-s010")
   sl = yt.SlicePlot(ds, "x", ("connect1", "diffused"))
   sl.zoom(0.75)
   sl.save()

Here, we plot the ``'diffused'`` variable, using a slice normal to the ``'x'`` direction,
through the meshed labelled by ``'connect1'``. By default, the slice goes through the
center of the domain. We have also zoomed out a bit to get a better view of the
resulting structure. To instead plot the ``'convected'`` variable, using a slice normal
to the ``'z'`` direction through the mesh labelled by ``'connect2'``, we do:

.. python-script::

   import yt

   ds = yt.load("MOOSE_sample_data/out.e-s010")
   sl = yt.SlicePlot(ds, "z", ("connect2", "convected"))
   sl.zoom(0.75)
   sl.save()

These slices are made by sampling the finite element solution at the points corresponding
to each pixel of the image. The ``'convected'`` and ``'diffused'`` variables are node-centered,
so this interpolation is performed by converting the sample point the reference coordinate
system of the element and evaluating the appropriate shape functions. You can also
plot element-centered fields:

.. python-script::

   import yt

   ds = yt.load("MOOSE_sample_data/out.e-s010")
   sl = yt.SlicePlot(ds, "y", ("connect1", "conv_indicator"))
   sl.zoom(0.75)
   sl.save()

We can also annotate the mesh lines, as follows:

.. python-script::

   import yt

   ds = yt.load("MOOSE_sample_data/out.e-s010")
   sl = yt.SlicePlot(ds, "z", ("connect1", "diffused"))
   sl.annotate_mesh_lines(plot_args={"color": "black"})
   sl.zoom(0.75)
   sl.save()

The ``plot_args`` parameter is a dictionary of keyword arguments that will be passed
to matplotlib. It can be used to control the mesh line color, thickness, etc...

The above examples all involve 8-node hexahedral mesh elements. Here is another example from
a dataset that uses 6-node wedge elements:

.. python-script::

   import yt

   ds = yt.load("MOOSE_sample_data/wedge_out.e")
   sl = yt.SlicePlot(ds, 2, ("connect2", "diffused"))
   sl.save()

Slices can also be used to examine 2D unstructured mesh datasets, but the
slices must be taken to be normal to the ``'z'`` axis, or you'll get an error. Here is
an example using another MOOSE dataset that uses triangular mesh elements:

.. python-script::

   import yt

   ds = yt.load("MOOSE_sample_data/out.e")
   sl = yt.SlicePlot(ds, 2, ("connect1", "nodal_aux"))
   sl.save()

You may run into situations where you have a variable you want to visualize that
exists on multiple mesh blocks. To view the variable on ``all`` mesh blocks,
simply pass ``all`` as the first argument of the field tuple:

.. python-script::

   import yt

   ds = yt.load("MultiRegion/two_region_example_out.e", step=-1)
   sl = yt.SlicePlot(ds, "z", ("all", "diffused"))
   sl.save()

.. _particle-plotting-workarounds:

Additional Notes for Plotting Particle Data
-------------------------------------------

Below are some important caveats to note when visualizing particle data.

1. Off axis slice plotting is not available for any particle data.
   However, axis-aligned slice plots (as described in :ref:`slice-plots`)
   will work.

2. Off axis projections (as in :ref:`off-axis-projection`) will only work
   for SPH particles, i.e., particles that have a defined smoothing length.

Two workaround methods are available for plotting non-SPH particles with off-axis
projections.

1. :ref:`smooth-non-sph` - this method involves extracting particle data to be
   reloaded with :ref:`~yt.loaders.load_particles` and using the
   :ref:`~yt.frontends.stream.data_structures.add_SPH_fields` function to
   create smoothing lengths. This works well for relatively small datasets,
   but is not parallelized and may take too long for larger data.

2. Plot from a saved
   :class:`~yt.data_objects.construction_data_containers.YTCoveringGrid`,
   :class:`~yt.data_objects.construction_data_containers.YTSmoothedCoveringGrid`,
   or :class:`~yt.data_objects.construction_data_containers.YTArbitraryGrid`
   dataset.

This second method is illustrated below. First, construct one of the grid data
objects listed above. Then, use the
:func:`~yt.data_objects.data_containers.YTDataContainer.save_as_dataset`
function (see :ref:`saving_data`) to save a deposited particle field
(see :ref:`deposited-particle-fields`) as a reloadable dataset. This dataset
can then be loaded and visualized using both off-axis projections and slices.
Note, the change in the field name from ``("deposit", "nbody_mass")`` to
``("grid", "nbody_mass")`` after reloading.

.. python-script::

   import yt

   ds = yt.load("gadget_cosmology_plus/snap_N128L16_132.hdf5")
   # create a 128^3 covering grid over the entire domain
   L = 7
   cg = ds.covering_grid(level=L, left_edge=ds.domain_left_edge, dims=[2**L]*3)

   fn = cg.save_as_dataset(fields=[("deposit", "nbody_mass")])
   ds_grid = yt.load(fn)
   p = yt.OffAxisProjectionPlot(ds_grid, [1, 1, 1], ("grid", "nbody_mass"))
   p.save()

Plot Customization: Recentering, Resizing, Colormaps, and More
--------------------------------------------------------------

You can customize each of the four plot types above in identical ways.  We'll go
over each of the customizations methods below.  For each of the examples below we
will modify the following plot.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.save()

Panning and zooming
~~~~~~~~~~~~~~~~~~~

There are three methods to dynamically pan around the data.

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.pan` accepts x and y
deltas.

.. python-script::

   import yt
   from yt.units import kpc

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.pan((2 * kpc, 2 * kpc))
   slc.save()

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.pan_rel` accepts deltas
in units relative to the field of view of the plot.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.pan_rel((0.1, -0.1))
   slc.save()

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.zoom` accepts a factor to zoom in by.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.zoom(2)
   slc.save()

Set axes units
~~~~~~~~~~~~~~

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_axes_unit` allows the customization of
the axes unit labels.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_axes_unit("Mpc")
   slc.save()

The same result could have been accomplished by explicitly setting the ``width``
to ``(.01, 'Mpc')``.

Set image units
~~~~~~~~~~~~~~~

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_axes_unit` allows
the customization of the units used for the image and colorbar.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_unit("density", "Msun/pc**3")
   slc.save()

If the unit you would like to convert to needs an equivalency, this can be
specified via the ``equivalency`` keyword argument of ``set_unit``. For
example, let's make a plot of the temperature field, but present it using
an energy unit instead of a temperature unit:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "temperature", width=(10, "kpc"))
   slc.set_unit("temperature", "keV", equivalency="thermal")
   slc.save()

Set the plot center
~~~~~~~~~~~~~~~~~~~

The :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_center`
function accepts a new center for the plot, in code units.  New centers must be
two element tuples.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_center((0.5, 0.503))
   slc.save()

Flipping the plot view axes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, all :class:`~yt.visualization.plot_window.PlotWindow` objects plot
with the assumption that the eastern direction on the plot forms a right handed
coordinate system with the ``normal`` and ``north_vector`` for the system, whether
explicitly or implicitly defined. This setting can be toggled or explicitly defined
by the user at initialization:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   # slicing with non right-handed coordinates
   slc = yt.SlicePlot(ds, "x", "velocity_x", right_handed=False)
   slc.annotate_title("Not Right Handed")
   slc.save("NotRightHanded.png")

   # switching to right-handed coordinates
   slc.toggle_right_handed()
   slc.annotate_title("Right Handed")
   slc.save("Standard.png")

.. _hiding-colorbar-and-axes:

Hiding the Colorbar and Axis Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~yt.visualization.plot_window.PlotWindow` class has functions
attached for hiding/showing the colorbar and axes.  This allows for making
minimal plots that focus on the data:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.hide_colorbar()
   slc.hide_axes()
   slc.save()

See the cookbook recipe :ref:`show-hide-axes-colorbar` and the full function
description :class:`~yt.visualization.plot_window.PlotWindow` for more
information.

Fonts
~~~~~

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_font` allows font
customization.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_font({"family": "sans-serif", "style": "italic", "weight": "bold", "size": 24})
   slc.save()

Colormaps
~~~~~~~~~

Each of these functions accept two arguments.  In all cases the first argument
is a field name.  This makes it possible to use different custom colormaps for
different fields tracked by the plot object.

To change the colormap for the plot, call the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_cmap` function.
Use any of the colormaps listed in the :ref:`colormaps` section.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_cmap("density", "RdBu_r")
   slc.save()

The :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_log` function
accepts a field name and a boolean.  If the boolean is ``True``, the colormap
for the field will be log scaled.  If it is ``False`` the colormap will be
linear.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_log("density", False)
   slc.save()

Specifically, a field containing both positive and negative values can be plotted
with symlog scale, by setting the boolean to be ``True`` and either providing an extra
parameter ``linthresh`` or setting ``symlog_auto = True``. In the region around zero
(when the log scale approaches to infinity), the linear scale will be applied to the
region ``(-linthresh, linthresh)`` and stretched relative to the logarithmic range.
In some cases, if yt detects zeros present in the dataset and the user has selected
``log`` scaling, yt automatically switches to ``symlog`` scaling and automatically
chooses a ``linthresh`` value to avoid errors.  This is the same behavior you can
achieve by setting the keyword ``symlog_auto`` to ``True``. In these cases, yt will
choose the smallest non-zero value in a dataset to be the ``linthresh`` value.
As an example,

.. python-script::

   import yt

   ds = yt.load_sample("FIRE_M12i_ref11")
   p = yt.ProjectionPlot(ds, "x", ("gas", "density"))
   p.set_log(("gas", "density"), True, symlog_auto=True)
   p.save()

Symlog is very versatile, and will work with positive or negative dataset ranges.
Here is an example using symlog scaling to plot a postive field with a linear range of
``(0, linthresh)``.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "x-velocity", width=(30, "kpc"))
   slc.set_log("x-velocity", True, linthresh=1.0e1)
   slc.save()

The :meth:`~yt.visualization.plot_container.ImagePlotContainer.set_background_color`
function accepts a field name and a color (optional). If color is given, the function
will set the plot's background color to that. If not, it will set it to the bottom
value of the color map.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(1.5, "Mpc"))
   slc.set_background_color("density")
   slc.save("bottom_colormap_background")
   slc.set_background_color("density", color="black")
   slc.save("black_background")

If you would like to change the background for a plot and also hide the axes,
you will need to make use of the ``draw_frame`` keyword argument for the ``hide_axes`` function. If you do not use this keyword argument, the call to
``set_background_color`` will have no effect. Here is an example illustrating how to use the ``draw_frame`` keyword argument for ``hide_axes``:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   field = ("deposit", "all_density")
   slc = yt.ProjectionPlot(ds, "z", field, width=(1.5, "Mpc"))
   slc.set_background_color(field)
   slc.hide_axes(draw_frame=True)
   slc.hide_colorbar()
   slc.save("just_image")

Lastly, the :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_zlim`
function makes it possible to set a custom colormap range.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_zlim("density", 1e-30, 1e-25)
   slc.save()

Annotations
~~~~~~~~~~~

A slice object can also add annotations like a title, an overlying
quiver plot, the location of grid boundaries, halo-finder annotations,
and many other annotations, including user-customizable annotations.
For example:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.annotate_grids()
   slc.save()

will plot the density field in a 10 kiloparsec slice through the
z-axis centered on the highest density point in the simulation domain.
Before saving the plot, the script annotates it with the grid
boundaries, which are drawn as lines in the plot, with colors going
from black to white depending on the AMR level of the grid.

Annotations are described in :ref:`callbacks`.

Set the size and resolution of the plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To set the size of the plot, use the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_figure_size` function.  The argument
is the size of the longest edge of the plot in inches.  View the full resolution
image to see the difference more clearly.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_figure_size(10)
   slc.save()

To change the resolution of the image, call the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_buff_size` function.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_buff_size(1600)
   slc.save()

Also see cookbook recipe :ref:`image-resolution-primer` for more information
about the parameters that determine the resolution of your images.

Turning off minorticks
~~~~~~~~~~~~~~~~~~~~~~

By default minorticks for the x and y axes are turned on.
The minorticks may be removed using the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_minorticks`
function, which either accepts a specific field name including the 'all' alias
and the desired state for the plot as 'on' or 'off'. There is also an analogous
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_colorbar_minorticks`
function for the colorbar axis.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, "z", "density", width=(10, "kpc"))
   slc.set_minorticks("all", False)
   slc.set_colorbar_minorticks("all", False)
   slc.save()


.. _matplotlib-customization:

Further customization via matplotlib
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each :class:`~yt.visualization.plot_window.PlotWindow` object is really a
container for plots - one plot for each field specified in the list of fields
supplied when the plot object is created. The individual plots can be
accessed via the ``plots`` dictionary attached to each
:class:`~yt.visualization.plot_window.PlotWindow` object:

.. code-block:: python

    slc = SlicePlot(ds, 2, ["density", "temperature"])
    dens_plot = slc.plots["density"]

In this example ``dens_plot`` is an instance of
:class:`~yt.visualization.plot_window.WindowPlotMPL`, an object that wraps the
matplotlib
`figure <https://matplotlib.org/stable/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure>`_
and `axes <https://matplotlib.org/stable/api/axes_api.html#matplotlib.axes.Axes>`_
objects.  We can access these matplotlib primitives via attributes of
``dens_plot``.

.. code-block:: python

    figure = dens_plot.figure
    axes = dens_plot.axes
    colorbar_axes = dens_plot.cax

These are the
`figure <https://matplotlib.org/stable/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure>`_
and `axes <https://matplotlib.org/stable/api/axes_api.html#matplotlib.axes.Axes>`_
objects that control the actual drawing of the plot.  Arbitrary plot
customizations are possible by manipulating these objects.  See
:ref:`matplotlib-primitives` for an example.

.. _how-to-make-1d-profiles:

1D Profile Plots
----------------

1D profiles are used to calculate the average or the sum of a given quantity
with respect to a second quantity.  Two common examples are the "average density
as a function of radius" or "the total mass within a given set of density bins."
When created, they default to the average: in fact, they default to the average
as weighted by the total cell mass.  However, this can be modified to take
either the total value or the average with respect to a different quantity.

Profiles operate on :ref:`data objects <data-objects>`; they will take the
entire data contained in a sphere, a prism, an extracted region and so on, and
they will calculate and use that as input to their calculation.  To make a 1D
profile plot, create a (:class:`~yt.visualization.profile_plotter.ProfilePlot`)
object, supplying the data object, the field for binning, and a list of fields
to be profiled.

.. python-script::

   import yt
   from yt.units import kpc

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_galaxy = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 10 * kpc, 3 * kpc)
   plot = yt.ProfilePlot(my_galaxy, "density", ["temperature"])
   plot.save()

This will create a :class:`~yt.data_objects.selection_data_containers.YTDisk`
centered at [0.5, 0.5, 0.5], with a normal vector of [0.0, 0.0, 1.0], radius of
10 kiloparsecs and height of 3 kiloparsecs and will then make a plot of the
mass-weighted average temperature as a function of density for all of the gas
contained in the cylinder.

We could also have made a profile considering only the gas in a sphere.
For instance:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = ds.sphere([0.5, 0.5, 0.5], (100, "kpc"))
   plot = yt.ProfilePlot(my_sphere, "temperature", ["mass"], weight_field=None)
   plot.save()

Note that because we have specified the weighting field to be ``None``, the
profile plot will display the accumulated cell mass as a function of temperature
rather than the average. Also note the use of a ``(value, unit)`` tuple. These
can be used interchangeably with units explicitly imported from ``yt.units`` when
creating yt plots.

We can also accumulate along the bin field of a ``ProfilePlot`` (the bin field
is the x-axis in a ``ProfilePlot``, in the last example the bin field is
``Temperature``) by setting the ``accumulation`` keyword argument to ``True``.
The following example uses ``weight_field = None`` and ``accumulation = True`` to
generate a plot of the enclosed mass in a sphere:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = ds.sphere([0.5, 0.5, 0.5], (100, "kpc"))
   plot = yt.ProfilePlot(
       my_sphere, "radius", ["mass"], weight_field=None, accumulation=True
   )
   plot.save()

You can also access the data generated by profiles directly, which can be
useful for overplotting average quantities on top of phase plots, or for
exporting and plotting multiple profiles simultaneously from a time series.
The ``profiles`` attribute contains a list of all profiles that have been
made.  For each item in the list, the x field data can be accessed with ``x``.
The profiled fields can be accessed from the dictionary ``field_data``.

.. code-block:: python

   plot = ProfilePlot(my_sphere, "temperature", ["mass"], weight_field=None)
   profile = plot.profiles[0]
   # print the bin field, in this case temperature
   print(profile.x)
   # print the profiled mass field
   print(profile["mass"])

Other options, such as the number of bins, are also configurable. See the
documentation for :class:`~yt.visualization.profile_plotter.ProfilePlot` for
more information.

Overplotting Multiple 1D Profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is often desirable to overplot multiple 1D profile to show evolution
with time.  This is supported with the ``from_profiles`` class method.
1D profiles are created with the :func:`~yt.data_objects.profiles.create_profile`
method and then given to the ProfilePlot object.

.. python-script::

   import yt

   # Create a time-series object.
   es = yt.load_simulation("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo")
   es.get_time_series(redshifts=[5, 4, 3, 2, 1, 0])

   # Lists to hold profiles, labels, and plot specifications.
   profiles = []
   labels = []

   # Loop over each dataset in the time-series.
   for ds in es:
       # Create a data container to hold the whole dataset.
       ad = ds.all_data()
       # Create a 1d profile of density vs. temperature.
       profiles.append(
           yt.create_profile(
               ad,
               ["temperature"],
               fields=["mass"],
               weight_field=None,
               accumulation=True,
           )
       )

       # Add labels
       labels.append("z = %.2f" % ds.current_redshift)

   # Create the profile plot from the list of profiles.
   plot = yt.ProfilePlot.from_profiles(profiles, labels=labels)

   # Save the image.
   plot.save()


Customizing axis limits
~~~~~~~~~~~~~~~~~~~~~~~

By default the x and y limits for ``ProfilePlot`` are determined using the
:class:`~yt.data_objects.derived_quantities.Extrema` derived quantity.  If you
want to create a plot with custom axis limits, you have two options.

First, you can create a custom profile object using
:func:`~yt.data_objects.profiles.create_profile`.
This function accepts a dictionary of ``(max, min)`` tuples keyed to field names.

.. python-script::

    import yt
    import yt.units as u

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    sp = ds.sphere("m", 10 * u.kpc)
    profiles = yt.create_profile(
        sp,
        "temperature",
        "density",
        weight_field=None,
        extrema={"temperature": (1e3, 1e7), "density": (1e-26, 1e-22)},
    )
    plot = yt.ProfilePlot.from_profiles(profiles)
    plot.save()

You can also make use of the
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_xlim` and
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_ylim` functions to
customize the axes limits of a plot that has already been created.  Note that
calling ``set_xlim`` is much slower than calling ``set_ylim``.  This is because
``set_xlim`` must recreate the profile object using the specified extrema.
Creating a profile directly via :func:`~yt.data_objects.profiles.create_profile`
might be significantly faster.
Note that since there is only one bin field, ``set_xlim``
does not accept a field name as the first argument.

.. python-script::

   import yt
   import yt.units as u

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   sp = ds.sphere("m", 10 * u.kpc)
   plot = yt.ProfilePlot(sp, "temperature", "density", weight_field=None)
   plot.set_xlim(1e3, 1e7)
   plot.set_ylim("density", 1e-26, 1e-22)
   plot.save()


Customizing Units
~~~~~~~~~~~~~~~~~

Units for both the x and y axis can be controlled via the
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_unit` method.
Adjusting the plot units does not require recreating the histogram, so adjusting
units will always be inexpensive, requiring only an in-place unit conversion.

In the following example we create a plot of the average density in solar
masses per cubic parsec as a function of radius in kiloparsecs.

.. python-script::

    import yt
    import yt.units as u

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    sp = ds.sphere("m", 10 * u.kpc)
    plot = yt.ProfilePlot(sp, "radius", "density", weight_field=None)
    plot.set_unit("density", "msun/pc**3")
    plot.set_unit("radius", "kpc")
    plot.save()

Linear and Logarithmic Scaling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The axis scaling can be manipulated via the
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_log` function.  This
function accepts a field name and a boolean.  If the boolean is ``True``, the
field is plotted in log scale.  If ``False``, the field is plotted in linear
scale.

In the following example we create a plot of the average x velocity as a
function of radius.  Since the x component of the velocity vector can be
negative, we set the scaling to be linear for this field.

.. python-script::

   import yt
   import yt.units as u

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   sp = ds.sphere("m", 10 * u.kpc)
   plot = yt.ProfilePlot(sp, "radius", "x-velocity", weight_field=None)
   plot.set_log("x-velocity", False)
   plot.save()

Setting axis labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The axis labels can be manipulated via the
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_ylabel` and
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_xlabel` functions.  The
:meth:`~yt.visualization.profile_plotter.ProfilePlot.set_ylabel` function accepts a field name
and a string with the desired label. The :meth:`~yt.visualization.profile_plotter.ProfilePlot.set_xlabel`
function just accepts the desired label and applies this to all of the plots.

In the following example we create a plot of the average x-velocity and density as a
function of radius. The xlabel is set to "Radius", for all plots, and the ylabel is set to
"velocity in x direction" for the x-velocity plot.

.. python-script::

   import yt

   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   ad = ds.all_data()
   plot = yt.ProfilePlot(ad, "density", ["temperature", "velocity_x"], weight_field=None)
   plot.set_xlabel("Radius")
   plot.set_ylabel("velocity_x", "velocity in x direction")
   plot.save()

Adding plot title
~~~~~~~~~~~~~~~~~

Plot title can be set via the
:meth:`~yt.visualization.profile_plotter.ProfilePlot.annotate_title` function.
It accepts a string argument which is the plot title and an optional ``field`` parameter which specifies
the field for which plot title should be added. ``field`` could be a string or a list of string.
If ``field`` is not passed, plot title will be added for the fields.

In the following example we create a plot and set the plot title.

.. python-script::

   import yt

   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   ad = ds.all_data()
   plot = yt.ProfilePlot(ad, "density", ["temperature"], weight_field=None)
   plot.annotate_title("Temperature vs Density Plot")
   plot.save()

Another example where we create plots from profile. By specifying the fields we can add plot title to a
specific plot.

.. python-script::

   import yt

   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   sphere = ds.sphere("max", (1.0, "Mpc"))
   profiles = []
   profiles.append(yt.create_profile(sphere, ["radius"], fields=["density"], n_bins=64))
   profiles.append(
       yt.create_profile(sphere, ["radius"], fields=["dark_matter_density"], n_bins=64)
   )
   plot = yt.ProfilePlot.from_profiles(profiles)
   plot.annotate_title("Plot Title: Density", "density")
   plot.annotate_title("Plot Title: Dark Matter Density", "dark_matter_density")
   plot.save()

Here, ``plot.annotate_title("Plot Title: Density", "density")`` will only set the plot title for the ``"density"``
field. Thus, allowing us the option to have different plot titles for different fields.


Annotating plot with text
~~~~~~~~~~~~~~~~~~~~~~~~~

Plots can be annotated at a desired (x,y) co-ordinate using :meth:`~yt.visualization.profile_plotter.ProfilePlot.annotate_text` function.
This function accepts the x-position, y-position, a text string to
be annotated in the plot area, and an optional list of fields for annotating plots with the specified field.
Furthermore, any keyword argument accepted by the matplotlib ``axes.text`` function could also be passed which will can be useful to change fontsize, text-alignment, text-color or other such properties of annotated text.

In the following example we create a plot and add a simple annotation.

.. python-script::

   import yt

   ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
   ad = ds.all_data()
   plot = yt.ProfilePlot(ad, "density", ["temperature"], weight_field=None)
   plot.annotate_text(1e-30, 1e7, "Annotated Text")
   plot.save()

To add annotations to a particular set of fields we need to pass in the list of fields as follows:

.. code-block:: python

   plot.annotate_text(1e-30, 1e7, "Annotation", ["field1", "field2"])


To change the text annotated text properties, we need to pass the matplotlib ``axes.text`` arguments as follows:

.. code-block:: python

  plot.annotate_text(
      1e-30,
      1e7,
      "Annotation",
      fontsize=20,
      bbox=dict(facecolor="red", alpha=0.5),
      horizontalalignment="center",
      verticalalignment="center",
  )

The above example will set the fontsize of annotation to 20, add a bounding box of red color and center align
horizontally and vertically. The is just an example to modify the text properties, for further options please check
`matplotlib.axes.Axes.text <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>`_.

Altering Line Properties
~~~~~~~~~~~~~~~~~~~~~~~~

Line properties for any and all of the profiles can be changed with the
:func:`~yt.visualization.profile_plotter.set_line_property` function.
The two arguments given are the line property and desired value.

.. code-block:: python

    plot.set_line_property("linestyle", "--")

With no additional arguments, all of the lines plotted will be altered.  To
change the property of a single line, give also the index of the profile.

.. code-block:: python

    # change only the first line
    plot.set_line_property("linestyle", "--", 0)

.. _how-to-1d-unstructured-mesh:

1D Line Sampling
----------------

YT has the ability to sample datasets along arbitrary lines
and plot the result. You must supply five arguments to the ``LinePlot``
class. They are enumerated below:

1. Dataset
2. A list of fields or a single field you wish to plot
3. The starting point of the sampling line. This should be an n-element list, tuple,
   ndarray, or YTArray with the elements corresponding to the coordinates of the
   starting point. (n should equal the dimension of the dataset)
4. The ending point of the sampling line. This should also be an n-element list, tuple,
   ndarray, or YTArray with the elements corresponding to the coordinates of the
   ending point.
5. The number of sampling points along the line, e.g. if 1000 is specified, then
   data will be sampled at 1000 points evenly spaced between the starting and
   ending points.

The below code snippet illustrates how this is done:

.. code-block:: python

   ds = yt.load("SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e", step=-1)
   plot = yt.LinePlot(ds, [("all", "v"), ("all", "u")], (0, 0, 0), (0, 1, 0), 1000)
   plot.save()

If working in a Jupyter Notebook, ``LinePlot`` also has the ``show()`` method.

You can add a legend to a 1D sampling plot. The legend process takes two steps:

1. When instantiating the ``LinePlot``, pass a dictionary of
   labels with keys corresponding to the field names
2. Call the ``LinePlot`` ``annotate_legend`` method

X- and Y- axis units can be set with ``set_x_unit`` and ``set_unit`` methods
respectively. The below code snippet combines all the features we've discussed:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   plot = yt.LinePlot(ds, "density", [0, 0, 0], [1, 1, 1], 512)
   plot.annotate_legend("density")
   plot.set_x_unit("cm")
   plot.set_unit("density", "kg/cm**3")
   plot.save()

If a list of fields is passed to ``LinePlot``, yt will create a number of
individual figures equal to the number of different dimensional
quantities. E.g. if ``LinePlot`` receives two fields with units of "length/time"
and a field with units of "temperature", two different figures will be created,
one with plots of the "length/time" fields and another with the plot of the
"temperature" field. It is only necessary to call ``annotate_legend``
for one field of a multi-field plot to produce a legend containing all the
labels passed in the initial construction of the ``LinePlot`` instance. Example:

.. python-script::

   import yt

   ds = yt.load("SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e", step=-1)
   plot = yt.LinePlot(
       ds,
       [("all", "v"), ("all", "u")],
       [0, 0, 0],
       [0, 1, 0],
       100,
       field_labels={("all", "u"): r"v$_x$", ("all", "v"): r"v$_y$"},
   )
   plot.annotate_legend(("all", "u"))
   plot.save()

``LinePlot`` is a bit different from yt ray objects which are data
containers. ``LinePlot`` is a plotting class that may use yt ray objects to
supply field plotting information. However, perhaps the most important
difference to highlight between rays and ``LinePlot`` is that rays return data
elements that intersect with the ray and make no guarantee about the spacing
between data elements. ``LinePlot`` sampling points are guaranteed to be evenly
spaced. In the case of cell data where multiple points fall within the same
cell, the ``LinePlot`` object will show the same field value for each sampling
point that falls within the same cell.

.. _how-to-make-2d-profiles:

2D Phase Plots
--------------

2D phase plots function in much the same was as 1D phase plots, but with a
:class:`~yt.visualization.profile_plotter.PhasePlot` object.  Much like 1D
profiles, 2D profiles (phase plots) are best thought of as plotting a
distribution of points, either taking the average or the accumulation in a bin.
The default behavior is to average, using the cell mass as the weighting,
but this behavior can be controlled through the ``weight_field`` parameter.
For example, to generate a 2D distribution of mass enclosed in density and
temperature bins, you can do:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = ds.sphere("c", (50, "kpc"))
   plot = yt.PhasePlot(
       my_sphere, "density", "temperature", ["mass"], weight_field=None
   )
   plot.save()

If you would rather see the average value of a field as a function of two other
fields, leave off the ``weight_field`` argument, and it will average by
the cell mass.  This would look
something like:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = ds.sphere("c", (50, "kpc"))
   plot = yt.PhasePlot(my_sphere, "density", "temperature", ["H_fraction"])
   plot.save()

Customizing Phase Plots
~~~~~~~~~~~~~~~~~~~~~~~

Similarly to 1D profile plots, :class:`~yt.visualization.profile_plotter.PhasePlot`
can be customized via ``set_unit``,
``set_xlim``, ``set_ylim``, and ``set_zlim``.  The following example illustrates
how to manipulate these functions. :class:`~yt.visualization.profile_plotter.PhasePlot`
can also be customized in a similar manner as
:class:`~yt.visualization.plot_window.SlicePlot`, such as with ``hide_colorbar``
and ``show_colorbar``.

.. python-script::

   import yt

   ds = yt.load("sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art")
   center = ds.arr([64.0, 64.0, 64.0], "code_length")
   rvir = ds.quan(1e-1, "Mpccm/h")
   sph = ds.sphere(center, rvir)

   plot = yt.PhasePlot(sph, "density", "temperature", "mass", weight_field=None)
   plot.set_unit("density", "Msun/pc**3")
   plot.set_unit("mass", "Msun")
   plot.set_xlim(1e-5, 1e1)
   plot.set_ylim(1, 1e7)
   plot.save()

It is also possible to construct a custom 2D profile object and then use the
:meth:`~yt.visualization.profile_plotter.PhasePlot.from_profile` function to
create a ``PhasePlot`` using the profile object.
This will sometimes be faster, especially if you need custom x and y axes
limits.  The following example illustrates this workflow:

.. python-script::

   import yt

   ds = yt.load("sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art")
   center = ds.arr([64.0, 64.0, 64.0], "code_length")
   rvir = ds.quan(1e-1, "Mpccm/h")
   sph = ds.sphere(center, rvir)
   units = dict(density="Msun/pc**3", cell_mass="Msun")
   extrema = dict(density=(1e-5, 1e1), temperature=(1, 1e7))

   profile = yt.create_profile(
       sph,
       ["density", "temperature"],
       n_bins=[128, 128],
       fields=["mass"],
       weight_field=None,
       units=units,
       extrema=extrema,
   )

   plot = yt.PhasePlot.from_profile(profile)

   plot.save()

Probability Distribution Functions and Accumulation
---------------------------------------------------

Both 1D and 2D profiles which show the total of amount of some field, such as
mass, in a bin (done by setting the ``weight_field`` keyword to ``None``) can be
turned into probability distribution functions (PDFs) by setting the
``fractional`` keyword to ``True``.  When set to ``True``, the value in each bin
is divided by the sum total from all bins.  These can be turned into cumulative
distribution functions (CDFs) by setting the ``accumulation`` keyword to
``True``.  This will make it so that the value in any bin N is the cumulative
sum of all bins from 0 to N.  The direction of the summation can be reversed by
setting ``accumulation`` to ``-True``.  For ``PhasePlot``, the accumulation can
be set independently for each axis by setting ``accumulation`` to a list of
``True``/ ``-True`` /``False`` values.

.. _particle-plots:

Particle Plots
--------------

Slice and projection plots both provide a callback for over-plotting particle
positions onto gas fields. However, sometimes you want to plot the particle
quantities by themselves, perhaps because the gas fields are not relevant to
your use case, or perhaps because your dataset doesn't contain any gas fields
in the first place. Additionally, you may want to plot your particles with a
third field, such as particle mass or age,  mapped to a colorbar.
:class:`~yt.visualization.particle_plots.ParticlePlot` provides a convenient
way to do this in yt.

The easiest way to make a :class:`~yt.visualization.particle_plots.ParticlePlot`
is to use the convenience routine. This has the syntax:

.. code-block:: python

   p = yt.ParticlePlot(ds, "particle_position_x", "particle_position_y")
   p.save()

Here, ``ds`` is a dataset we've previously opened. The commands create a particle
plot that shows the x and y positions of all the particles in ``ds`` and save the
result to a file on the disk. The type of plot returned depends on the fields you
pass in; in this case, ``p`` will be an :class:`~yt.visualization.particle_plots.ParticleProjectionPlot`,
because the fields are aligned to the coordinate system of the simulation.
The above example is equivalent to the following:

.. code-block:: python

   p = yt.ParticleProjectionPlot(ds, "z")
   p.save()

Most of the callbacks the work for slice and projection plots also work for
:class:`~yt.visualization.particle_plots.ParticleProjectionPlot`.
For instance, we can zoom in:

.. code-block:: python

   p = yt.ParticlePlot(ds, "particle_position_x", "particle_position_y")
   p.zoom(10)
   p.save("zoom")

change the width:

.. code-block:: python

   p.set_width((500, "kpc"))

or change the axis units:

.. code-block:: python

   p.set_unit("particle_position_x", "Mpc")

Here is a full example that shows the simplest way to use
:class:`~yt.visualization.particle_plots.ParticlePlot`:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ParticlePlot(ds, "particle_position_x", "particle_position_y")
   p.save()

In the above examples, we are simply splatting particle x and y positions onto
a plot using some color. Colors can be applied to the plotted particles by
providing a ``z_field``, which will be summed along the line of sight in a manner
similar to a projection.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ParticlePlot(ds, "particle_position_x", "particle_position_y", "particle_mass")
   p.set_unit("particle_mass", "Msun")
   p.zoom(32)
   p.save()

Additionally, a ``weight_field`` can be given such that the value in each
pixel is the weighted average along the line of sight.

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ParticlePlot(
       ds,
       "particle_position_x",
       "particle_position_y",
       "particle_mass",
       weight_field="particle_ones",
   )
   p.set_unit("particle_mass", "Msun")
   p.zoom(32)
   p.save()

Note the difference in the above two plots. The first shows the
total mass along the line of sight. The density is higher in the
inner regions, and hence there are more particles and more mass along
the line of sight. The second plot shows the average mass per particle
along the line of sight. The inner region is dominated by low mass
star particles, whereas the outer region is comprised of higher mass
dark matter particles.

Both :class:`~yt.visualization.particle_plots.ParticleProjectionPlot` and
:class:`~yt.visualization.particle_plots.ParticlePhasePlot` objects
accept a ``deposition`` argument which controls the order of the "splatting"
of the particles onto the pixels in the plot. The default option, ``"ngp"``,
corresponds to the "Nearest-Grid-Point" (0th-order) method, which simply
finds the pixel the particle is located in and deposits 100% of the particle
or its plotted quantity into that pixel. The other option, ``"cic"``,
corresponds to the "Cloud-In-Cell" (1st-order) method, which linearly
interpolates the particle or its plotted quantity into the four nearest
pixels in the plot.

Here is a complete example that uses the ``particle_mass`` field
to set the colorbar and shows off some of the modification functions for
:class:`~yt.visualization.particle_plots.ParticleProjectionPlot`:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ParticlePlot(
       ds,
       "particle_position_x",
       "particle_position_y",
       "particle_mass",
       width=(0.5, 0.5),
   )
   p.set_unit("particle_mass", "Msun")
   p.zoom(32)
   p.annotate_title("Zoomed-in Particle Plot")
   p.save()

If the fields passed in to :class:`~yt.visualization.particle_plots.ParticlePlot`
do not correspond to a valid :class:`~yt.visualization.particle_plots.ParticleProjectionPlot`,
a :class:`~yt.visualization.particle_plots.ParticlePhasePlot` will be returned instead.
:class:`~yt.visualization.particle_plots.ParticlePhasePlot` is used to plot arbitrary particle
fields against each other, and do not support some of the callbacks available in
:class:`~yt.visualization.particle_plots.ParticleProjectionPlot` -
for instance, :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.pan` and
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.zoom` don't make much sense when of your axes is a position
and the other is a velocity. The modification functions defined for :class:`~yt.visualization.profile_plotter.PhasePlot`
should all work, however.

Here is an example of making a :class:`~yt.visualization.particle_plots.ParticlePhasePlot`
of ``particle_position_x`` versus ``particle_velocity_z``, with the ``particle_mass`` on the colorbar:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ParticlePlot(ds, "particle_position_x", "particle_velocity_z", ["particle_mass"])
   p.set_unit("particle_position_x", "Mpc")
   p.set_unit("particle_velocity_z", "km/s")
   p.set_unit("particle_mass", "Msun")
   p.save()

and here is one with the particle x and y velocities on the plot axes:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ParticlePlot(ds, "particle_velocity_x", "particle_velocity_y", "particle_mass")
   p.set_unit("particle_velocity_x", "km/s")
   p.set_unit("particle_velocity_y", "km/s")
   p.set_unit("particle_mass", "Msun")
   p.set_ylim(-400, 400)
   p.set_xlim(-400, 400)
   p.save()

If you want more control over the details of the :class:`~yt.visualization.particle_plots.ParticleProjectionPlot` or
:class:`~yt.visualization.particle_plots.ParticlePhasePlot`, you can always use these classes directly. For instance,
here is an example of using the ``depth`` argument to :class:`~yt.visualization.particle_plots.ParticleProjectionPlot`
to only plot the particles that live in a thin slice around the center of the
domain:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   p = yt.ParticleProjectionPlot(ds, 2, ["particle_mass"], width=(0.5, 0.5), depth=0.01)
   p.set_unit("particle_mass", "Msun")
   p.save()

and here is an example of using the ``data_source`` argument to :class:`~yt.visualization.particle_plots.ParticlePhasePlot`
to only consider the particles that lie within a 50 kpc sphere around the domain center:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   my_sphere = ds.sphere("c", (50.0, "kpc"))

   p = yt.ParticlePhasePlot(
       my_sphere, "particle_velocity_x", "particle_velocity_y", "particle_mass"
   )
   p.set_unit("particle_velocity_x", "km/s")
   p.set_unit("particle_velocity_y", "km/s")
   p.set_unit("particle_mass", "Msun")
   p.set_ylim(-400, 400)
   p.set_xlim(-400, 400)

   p.save()

:class:`~yt.visualization.particle_plots.ParticleProjectionPlot` objects also admit a ``density``
flag, which allows one to plot the surface density of a projected quantity. This simply divides
the quantity in each pixel of the plot by the area of that pixel. It also changes the label on the
colorbar to reflect the new units and the fact that it is a density. This may make most sense in
the case of plotting the projected particle mass, in which case you can plot the projected particle
mass density:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   p = yt.ParticleProjectionPlot(ds, 2, ["particle_mass"], width=(0.5, 0.5), density=True)
   p.set_unit("particle_mass", "Msun/kpc**2") # Note that the dimensions reflect the density flag
   p.save()

Finally, with 1D and 2D Profiles, you can create a :class:`~yt.data_objects.profiles.ParticleProfile`
object separately using the :func:`~yt.data_objects.profiles.create_profile` function, and then use it
create a :class:`~yt.visualization.particle_plots.ParticlePhasePlot` object using the
:meth:`~yt.visualization.particle_plots.ParticlePhasePlot.from_profile` method. In this example,
we have also used the ``weight_field`` argument to compute the average ``particle_mass`` in each
pixel, instead of the total:

.. python-script::

   import yt

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

   ad = ds.all_data()

   profile = yt.create_profile(
       ad,
       ["particle_velocity_x", "particle_velocity_y"],
       ["particle_mass"],
       n_bins=800,
       weight_field="particle_ones",
   )

   p = yt.ParticlePhasePlot.from_profile(profile)
   p.set_unit("particle_velocity_x", "km/s")
   p.set_unit("particle_velocity_y", "km/s")
   p.set_unit("particle_mass", "Msun")
   p.set_ylim(-400, 400)
   p.set_xlim(-400, 400)
   p.save()

Under the hood, the :class:`~yt.data_objects.profiles.ParticleProfile` class works a lot like a
:class:`~yt.data_objects.profiles.Profile2D` object, except that instead of just binning the
particle field, you can also use higher-order deposition functions like the cloud-in-cell
interpolant to spread out the particle quantities over a few cells in the profile. The
:func:`~yt.data_objects.profiles.create_profile` will automatically detect when all the fields
you pass in are particle fields, and return a :class:`~yt.data_objects.profiles.ParticleProfile`
if that is the case. For a complete description of the :class:`~yt.data_objects.profiles.ParticleProfile`
class please consult the reference documentation.

.. _interactive-plotting:

Interactive Plotting
--------------------

The best way to interactively plot data is through the IPython notebook.  Many
detailed tutorials on using the IPython notebook can be found at
:ref:`notebook-tutorial`. The simplest way to launch the notebook it is to
type:

.. code-block:: bash

   yt notebook

at the command line.  This will prompt you for a password (so that if you're on
a shared user machine no one else can pretend to be you!) and then spawn an
IPython notebook you can connect to.

If you want to see yt plots inline inside your notebook, you need only create a
plot and then call ``.show()`` and the image will appear inline:

.. notebook-cell::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = yt.ProjectionPlot(ds, "z", "density", center='m', width=(10,'kpc'),
                      weight_field='density')
   p.set_figure_size(5)
   p.show()

.. _saving_plots:

Saving Plots
------------

If you want to save your yt plots, you have a couple of options for customizing
the plot filenames. If you don't care what the filenames are, just calling the
``save`` method with no additional arguments usually suffices:

.. code-block:: python

   import yt

   ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
   slc = yt.SlicePlot(ds, "z", ["kT", "density"], width=(500.0, "kpc"))
   slc.save()

which will yield PNG plots with the filenames

.. code-block:: bash

   $ ls \*.png
   sloshing_nomag2_hdf5_plt_cnt_0100_Slice_z_density.png
   sloshing_nomag2_hdf5_plt_cnt_0100_Slice_z_kT.png

which has a general form of

.. code-block:: bash

   [dataset name]_[plot type]_[axis]_[field name].[suffix]

Calling ``save`` with a single argument or the ``name`` keyword argument
specifies an alternative name for the plot:

.. code-block:: python

   slc.save("bananas")

or

.. code-block:: python

   slc.save(name="bananas")

yields

.. code-block:: bash

   $ ls \*.png
   bananas_Slice_z_kT.png
   bananas_Slice_z_density.png

If you call ``save`` with a full filename with a file suffix, the plot
will be saved with that filename:

.. code-block:: python

   slc.save("sloshing.png")

since this will take any field and plot it with this filename, it is
typically only useful if you are plotting one field. If you want to
simply change the image format of the plotted file, use the ``suffix``
keyword:

.. code-block:: python

   slc.save(name="bananas", suffix="eps")

yielding

.. code-block:: bash

   $ ls *.eps
   bananas_Slice_z_kT.eps
   bananas_Slice_z_density.eps

.. _remaking-plots:

Remaking Figures from Plot Datasets
-----------------------------------

When working with datasets that are too large to be stored locally,
making figures just right can be cumbersome as it requires continuously
moving images somewhere they can be viewed.  However, image creation is
actually a two step process of first creating the projection, slice,
or profile object, and then converting that object into an actual image.
Fortunately, the hard part (creating slices, projections, profiles) can
be separated from the easy part (generating images).  The intermediate
slice, projection, and profile objects can be saved as reloadable
datasets, then handed back to the plotting machinery discussed here.

For slices and projections, the savable object is associated with the
plot object as ``data_source``.  This can be saved with the
:func:`~yt.data_objects.data_containers.save_as_dataset`` function.  For
more information, see :ref:`saving_data`.

.. code-block:: python

   p = yt.ProjectionPlot(ds, "x", "density", weight_field="density")
   fn = p.data_source.save_as_dataset()

This function will optionally take a ``filename`` keyword that follows
the same logic as discussed above in :ref:`saving_plots`.  The filename
to which the dataset was written will be returned.

Once saved, this file can be reloaded completely independently of the
original dataset and given back to the plot function with the same
arguments.  One can now continue to tweak the figure to one's liking.

.. code-block:: python

   new_ds = yt.load(fn)
   new_p = yt.ProjectionPlot(new_ds, "x", "density", weight_field="density")
   new_p.save()

The same functionality is available for profile and phase plots.  In
each case, a special data container, ``data``, is given to the plotting
functions.

For ``ProfilePlot``:

.. code-block:: python

   ad = ds.all_data()
   p1 = yt.ProfilePlot(ad, "density", "temperature", weight_field="mass")

   # note that ProfilePlots can hold a list of profiles
   fn = p1.profiles[0].save_as_dataset()

   new_ds = yt.load(fn)
   p2 = yt.ProfilePlot(new_ds.data, "density", "temperature", weight_field="mass")
   p2.save()

For ``PhasePlot``:

.. code-block:: python

   ad = ds.all_data()
   p1 = yt.PhasePlot(ad, "density", "temperature", "mass", weight_field=None)
   fn = p1.profile.save_as_dataset()

   new_ds = yt.load(fn)
   p2 = yt.PhasePlot(new_ds.data, "density", "temperature", "mass", weight_field=None)
   p2.save()

.. _eps-writer:

Publication-ready Figures
-------------------------

While the routines above give a convenient method to inspect and
visualize your data, publishers often require figures to be in PDF or
EPS format.  While the matplotlib supports vector graphics and image
compression in PDF formats, it does not support compression in EPS
formats.  The :class:`~yt.visualization.eps_writer.DualEPS` module
provides an interface with the `PyX <https://pyx-project.org/>`_,
which is a Python abstraction of the PostScript drawing model with a
LaTeX interface.  It is optimal for publications to provide figures
with vector graphics to avoid rasterization of the lines and text,
along with compression to produce figures that do not have a large
filesize.

.. note::
   PyX must be installed, which can be accomplished either manually
   with ``pip install pyx`` or with the install script by setting
   ``INST_PYX=1``. If you are using python2, you must install pyx
   version 0.12.1 with ``pip install pyx==0.12.1``, since that is
   the last version with python2 support.

This module can take any of the plots mentioned above and create an
EPS or PDF figure.  For example,

.. code-block:: python

    import yt.visualization.eps_writer as eps

    slc = yt.SlicePlot(ds, "z", "density")
    slc.set_width(25, "kpc")
    eps_fig = eps.single_plot(slc)
    eps_fig.save_fig("zoom", format="eps")
    eps_fig.save_fig("zoom-pdf", format="pdf")

The ``eps_fig`` object exposes all of the low-level functionality of
``PyX`` for further customization (see the `PyX documentation
<https://pyx-project.org/manual/>`_).  There are a few
convenience routines in ``eps_writer``, such as drawing a circle,

.. code-block:: python

    eps_fig.circle(radius=0.2, loc=(0.5, 0.5))
    eps_fig.sav_fig("zoom-circle", format="eps")

with a radius of 0.2 at a center of (0.5, 0.5), both of which are in
units of the figure's field of view.  The
:func:`~yt.visualization.eps_writer.multiplot_yt` routine also
provides a convenient method to produce multi-panel figures
from a PlotWindow.  For example,

.. code-block:: python

    import yt
    import yt.visualization.eps_writer as eps

    slc = yt.SlicePlot(
        ds, "z", ["density", "temperature", "pressure", "velocity_magnitude"]
    )
    slc.set_width(25, "kpc")
    eps_fig = eps.multiplot_yt(2, 2, slc, bare_axes=True)
    eps_fig.scale_line(0.2, "5 kpc")
    eps_fig.save_fig("multi", format="eps")

will produce a 2x2 panel figure with a scale bar indicating 5 kpc.
The routine will try its best to place the colorbars in the optimal
margin, but it can be overridden by providing the keyword
``cb_location`` with a dict of either ``right, left, top, bottom``
with the fields as the keys.

You can also combine slices, projections, and phase plots. Here is
an example that includes slices and phase plots:

.. code-block:: python

    from yt import PhasePlot, SlicePlot
    from yt.visualization.eps_writer import multiplot_yt

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    p1 = SlicePlot(ds, 0, "density")
    p1.set_width(10, "kpc")

    p2 = SlicePlot(ds, 0, "temperature")
    p2.set_width(10, "kpc")
    p2.set_cmap("temperature", "hot")

    sph = ds.sphere(ds.domain_center, (10, "kpc"))
    p3 = PhasePlot(sph, "radius", "density", "temperature", weight_field="mass")

    p4 = PhasePlot(sph, "radius", "density", "pressure", "mass")

    mp = multiplot_yt(
        2,
        2,
        [p1, p2, p3, p4],
        savefig="yt",
        shrink_cb=0.9,
        bare_axes=False,
        yt_nocbar=False,
        margins=(0.5, 0.5),
    )

    mp.save_fig("multi_slice_phase")


.. _how-to-make-plots:

How to Make Plots
=================

In this section we explain how to use yt to create visualizations
of simulation data, derived fields, and the data produced by yt
analysis objects.  For details about the data extraction and
algorithms used to produce the image and analysis data, please see the
yt `method paper
<http://adsabs.harvard.edu/abs/2011ApJS..192....9T>`_.  There are also
many example scripts in :ref:`cookbook`.

The :class:`~yt.visualization.plot_window.PlotWindow` interface is useful for
taking a quick look at simulation outputs.  Simple mechanisms exist for making 
plots of slices, projections, 1D profiles, and 2D profiles (phase plots), all of 
which are described below.

.. _simple-inspection:

Slices & Projections
--------------------

If you need to take a quick look at a single simulation output, yt
provides the :class:`~yt.visualization.plot_window.PlotWindow` interface for 
generating annotated 2D visualizations of simulation data.  You can create a 
:class:`~yt.visualization.plot_window.PlotWindow` plot by
supplying a dataset, a list of fields to plot, and a plot center to
create a :class:`~yt.visualization.plot_window.AxisAlignedSlicePlot`,
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
`figure <http://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure>`
and `axes <http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes>` 
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

    slc = yt.SlicePlot(ds, 'z', 'density')
    slc.save()

These two commands will create a slice object and store it in a variable we've
called ``slc``.  Since this plot is aligned with the simulation coordinate
system, ``slc`` is an instance of
:class:`~yt.visualization.plot_window.AxisAlignedSlicePlot`. We then call the
``save()`` function, which automatically saves the plot in png image format with
an automatically generated filename.  If you don't want the slice object to
stick around, you can accomplish the same thing in one line:

.. code-block:: python
   
    yt.SlicePlot(ds, 'z', 'density').save()

It's nice to keep the slice object around if you want to modify the plot.  By
default, the plot width will be set to the size of the simulation box.  To zoom
in by a factor of ten, you can call the zoom function attached to the slice
object:

.. code-block:: python

    slc = yt.SlicePlot(ds, 'z', 'density')
    slc.zoom(10)
    slc.save('zoom')

This will save a new plot to disk with a different filename - prepended with
'zoom' instead of the name of the dataset. If you want to set the width
manually, you can do that as well. For example, the following sequence of
commands will create a slice, set the width of the plot to 10 kiloparsecs, and
save it to disk.

.. code-block:: python

    from yt.units import kpc
    slc = yt.SlicePlot(ds, 'z', 'density')
    slc.set_width(10*kpc)
    slc.save('10kpc')

The plot width can be specified independently along the x and y direction by
passing a tuple of widths.  An individual width can also be represented using a
``(value, unit)`` tuple.  The following sequence of commands all equivalently
set the width of the plot to 200 kiloparsecs in the ``x`` and ``y`` direction.

.. code-block:: python

    from yt.units import kpc
    slc.set_width(200*kpc)
    slc.set_width((200, 'kpc'))
    slc.set_width((200*kpc, 200*kpc))

The ``SlicePlot`` also optionally accepts the coordinate to center the plot on
and the width of the plot:

.. code-block:: python

    yt.SlicePlot(ds, 'z', 'density', center=[0.2, 0.3, 0.8],
                 width = (10,'kpc')).save()

Note that, by default,
:class:`~yt.visualization.plot_window.SlicePlot` shifts the
coordinates on the axes such that the origin is at the center of the
slice.  To instead use the coordinates as defined in the dataset, use
the optional argument: ``origin="native"``

If supplied without units, the center is assumed by in code units.  There are also
the following alternative options for the `center` keyword:

* ``"center"``, ``"c"``: the domain center
* ``"max"``, ``"m"``: the position of the maximum density
* ``("min", field)``: the position of the minimum of ``field``
* ``("max", field)``: the position of the maximum of ``field``

where for the last two objects any spatial field, such as ``"density"``,
``"velocity_z"``,
etc., may be used, e.g. ``center=("min","temperature")``.

Here is an example that combines all of the options we just discussed.

.. python-script::

   import yt
   from yt.units import kpc
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', center=[0.5, 0.5, 0.5],
                      width=(20,'kpc'))
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
   slc = yt.SlicePlot(ds, 'z', 'pressure', center='c')
   slc.save()
   slc.zoom(30)
   slc.save('zoom')

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

.. _off-axis-slices:

Off Axis Slices
~~~~~~~~~~~~~~~

Off axis slice plots can be generated in much the same way as
grid-aligned slices.  Off axis slices use
:class:`~yt.data_objects.selection_data_containers.YTCuttingPlaneBase` to slice
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
   L = [1,1,0] # vector normal to cutting plane
   north_vector = [-1,1,0]
   cut = yt.SlicePlot(ds, L, 'density', width=(25, 'kpc'),
                      north_vector=north_vector)
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
   prj = yt.ProjectionPlot(ds, 2, 'temperature', width=25*kpc,
                           weight_field='density')
   prj.save()

will create a density-weighted projection of the temperature field along the x
axis, plot it, and then save the plot to a png image file.

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

Internally, off axis projections are created using :ref:`the-camera-interface`
by applying the
:class:`~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction`.
In this use case, the volume renderer casts a set of plane parallel rays, one
for each pixel in the image.  The data values along each ray are summed,
creating the final image buffer.

.. _off-axis-projection-function:

To avoid manually creating a camera and setting the transfer
function, yt provides the :func:`~yt.visualization.volume_rendering.camera.off_axis_projection`
function, which wraps the camera interface to create an off axis
projection image buffer.  These images can be saved to disk or
used in custom plots.  This snippet creates an off axis
projection through a simulation.

.. python-script::

   import yt
   import numpy as np
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1,1,0] # vector normal to cutting plane
   north_vector = [-1,1,0]
   W = [0.02, 0.02, 0.02]
   c = [0.5, 0.5, 0.5]
   N = 512
   image = yt.off_axis_projection(ds, c, L, W, N, "density")
   yt.write_image(np.log10(image), "%s_offaxis_projection.png" % ds)

Here, ``W`` is the width of the projection in the x, y, *and* z
directions.

One can also generate generate annotated off axis projections
using
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`. These
plots can be created in much the same way as an
``OffAxisSlicePlot``, requiring only an open dataset, a direction
to project along, and a field to project.  For example:

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1,1,0] # vector normal to cutting plane
   north_vector = [-1,1,0]
   prj = yt.OffAxisProjectionPlot(ds,L,'density',width=(25, 'kpc'),
                                  north_vector=north_vector)
   prj.save()

OffAxisProjectionPlots can also be created with a number of
keyword arguments, as described in
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`

Plot Customization: Recentering, Resizing, Colormaps, and More
--------------------------------------------------------------

You can customize each of the four plot types above in identical ways.  We'll go
over each of the customizations methods below.  For each of the examples below we
will modify the following plot.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
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
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.pan((2*kpc, 2*kpc))
   slc.save()

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.pan_rel` accepts deltas 
in units relative to the field of view of the plot.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.pan_rel((0.1, -0.1))
   slc.save()

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.zoom` accepts a factor to zoom in by.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.zoom(2)
   slc.save()

Set axes units
~~~~~~~~~~~~~~

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_axes_unit` allows the customization of
the axes unit labels.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_axes_unit('Mpc')
   slc.save()

The same result could have been accomplished by explicitly setting the ``width``
to ``(.01, 'Mpc')``.

Set the plot center
~~~~~~~~~~~~~~~~~~~

The :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_center`
function accepts a new center for the plot, in code units.  New centers must be
two element tuples.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_center((0.5, 0.503))
   slc.save()


.. _hiding-colorbar-and-axes:

Hiding the Colorbar and Axis Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~yt.visualization.plot_window.PlotWindow` class has functions
attached for hiding/showing the colorbar and axes.  This allows for making
minimal plots that focus on the data:

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.hide_colorbar()
   slc.hide_axes()
   slc.save()

See the cookbook recipe :ref:`show-hide-axes-colorbar` and the 
`full function description ~yt.visualization.plot_window.PlotWindow` for more
information.

Fonts
~~~~~

:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_font` allows font
costomization.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_font({'family': 'sans-serif', 'style': 'italic',
                 'weight': 'bold', 'size': 24})
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
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_cmap('density', 'RdBu_r')
   slc.save()

The :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_log` function
accepts a field name and a boolean.  If the boolean is ``True``, the colormap
for the field will be log scaled.  If it is ``False`` the colormap will be
linear.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_log('density', False)
   slc.save()

Specifically, a field containing both positive and negative values can be plotted
with symlog scale, by seting the boolean to be ``True`` and providing an extra
parameter ``linthresh``. In the region around zero (when the log scale approaches
to infinity), the linear scale will be applied to the region ``(-linthresh, linthresh)``
and stretched relative to the logarithmic range. You can also plot a positive field 
under symlog scale with the linear range of ``(0, linthresh)``.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'x-velocity', width=(30,'kpc'))
   slc.set_log('x-velocity', True, linthresh=1.e1)
   slc.save()

Lastly, the :meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_zlim`
function makes it possible to set a custom colormap range.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_zlim('density', 1e-30, 1e-25)
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
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.annotate_grids()
   slc.save()

will plot the density field in a 10 kiloparsec slice through the
z-axis centered on the highest density point in the simulation domain.
Before saving the plot, the script annotates it with the grid
boundaries, which are drawn as lines in the plot, with colors going
from black to white depending on the AMR level of the grid.

Annotations are described in :ref:`callbacks`.

Set the size of the plot
~~~~~~~~~~~~~~~~~~~~~~~~

To set the size of the plot, use the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_figure_size` function.  The argument
is the size of the longest edge of the plot in inches.  View the full resolution
image to see the difference more clearly.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_figure_size(10)
   slc.save()

To change the resolution of the image, call the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_buff_size` function.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_buff_size(1600)
   slc.save()

Turning off minorticks
~~~~~~~~~~~~~~~~~~~~~~

By default minorticks for the x and y axes are turned on.
The minorticks may be removed using the
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_minorticks`
function, which either accepts a specific field name including the 'all' alias
and the desired state for the plot as 'on' or 'off'. There is also an analogous
:meth:`~yt.visualization.plot_window.AxisAlignedSlicePlot.set_cbar_minorticks`
function for the colorbar axis.

.. python-script::

   import yt
   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = yt.SlicePlot(ds, 'z', 'density', width=(10,'kpc'))
   slc.set_minorticks('all', 'off')
   slc.set_cbar_minorticks('all', 'off')
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

    slc = SlicePlot(ds, 2, ['density', 'temperature']
    dens_plot = slc.plots['density']

In this example ``dens_plot`` is an instance of
:class:`~yt.visualization.plot_window.WindowPlotMPL`, an object that wraps the
matplotlib 
`figure <http://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure>` 
and `axes <http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes>` 
objects.  We can access these matplotlib primitives via attributes of 
``dens_plot``.  

.. code-block:: python

    figure = dens_plot.figure
    axes = dens_plot.axes
    colorbar_axes = dens_plot.cax

These are the 
`figure <http://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure>`, 
and `axes <http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes>` 
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
   my_galaxy = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 10*kpc, 3*kpc)
   plot = yt.ProfilePlot(my_galaxy, "density", ["temperature"])
   plot.save()

This will create a :class:`~yt.data_objects.selection_data_containers.YTDiskBase`
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
   plot = yt.ProfilePlot(my_sphere, "temperature", ["cell_mass"],
                         weight_field=None)
   plot.save()

Note that because we have specified the weighting field to be ``None``, the
profile plot will display the accumulated cell mass as a function of temperature
rather than the average. Also note the use of a ``(value, unit)`` tuple. These
can be used interchangably with units explicitly imported from ``yt.units`` when
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
   plot = yt.ProfilePlot(my_sphere, "radius", ["cell_mass"],
                         weight_field=None, accumulation=True)
   plot.save()

You can also access the data generated by profiles directly, which can be
useful for overplotting average quantities on top of phase plots, or for
exporting and plotting multiple profiles simultaneously from a time series.
The ``profiles`` attribute contains a list of all profiles that have been 
made.  For each item in the list, the x field data can be accessed with ``x``.  
The profiled fields can be accessed from the dictionary ``field_data``.

.. code-block:: python

   plot = ProfilePlot(my_sphere, "temperature", ["cell_mass"],
                      weight_field=None)
   profile = plot.profiles[0]
   # print the bin field, in this case temperature
   print profile.x
   # print the profiled cell_mass field
   print profile['cell_mass']

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
   es = yt.simulation("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo")
   es.get_time_series(redshifts=[5, 4, 3, 2, 1, 0])


   # Lists to hold profiles, labels, and plot specifications.
   profiles = []
   labels = []

   # Loop over each dataset in the time-series.
   for ds in es:
       # Create a data container to hold the whole dataset.
       ad = ds.all_data()
       # Create a 1d profile of density vs. temperature.
       profiles.append(yt.create_profile(ad, ["temperature"], 
                                         fields=["cell_mass"],
                                         weight_field=None,
                                         accumulation=True))
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
    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    sp = ds.sphere('m', 10*u.kpc)
    profiles = yt.create_profile(sp, "temperature", "density",
                                 weight_field=None, 
                                 extrema={'temperature': (1e3, 1e7),
                                          'density': (1e-26, 1e-22)})
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
   ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
   sp = ds.sphere('m', 10*u.kpc)
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
    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    sp = ds.sphere('m', 10*u.kpc)
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
   ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
   sp = ds.sphere('m', 10*u.kpc)
   plot = yt.ProfilePlot(sp, "radius", "x-velocity", weight_field=None)
   plot.set_log("x-velocity", False)
   plot.save()

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
   plot = yt.PhasePlot(my_sphere, "density", "temperature", ["cell_mass"],
                       weight_field=None)
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
how to manipulate these functions.

.. python-script::

   import yt
   ds = yt.load("sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art")
   center = ds.arr([64.0, 64.0, 64.0], 'code_length')
   rvir = ds.quan(1e-1, "Mpccm/h")
   sph = ds.sphere(center, rvir)

   plot = yt.PhasePlot(sph, "density", "temperature", "cell_mass",
                       weight_field=None)
   plot.set_unit('density', 'Msun/pc**3')
   plot.set_unit('cell_mass', 'Msun')
   plot.set_xlim(1e-5,1e1)
   plot.set_ylim(1,1e7)
   plot.save()

It is also possible to construct a custom 2D profile object and then use the
:meth:`~yt.visualization.profile_plotter.PhasePlot.from_profile` function to 
create a ``PhasePlot`` using the profile object.
This will sometimes be faster, especially if you need custom x and y axes
limits.  The following example illustrates this workflow:

.. python-script::

   import yt
   ds = yt.load("sizmbhloz-clref04SNth-rs9_a0.9011/sizmbhloz-clref04SNth-rs9_a0.9011.art")
   center = ds.arr([64.0, 64.0, 64.0], 'code_length')
   rvir = ds.quan(1e-1, "Mpccm/h")
   sph = ds.sphere(center, rvir)
   units = dict(density='Msun/pc**3', cell_mass='Msun')
   extrema = dict(density=(1e-5, 1e1), temperature=(1, 1e7))

   profile = yt.create_profile(sph, ['density', 'temperature'],
                               n_bins=[128, 128], fields=['cell_mass'],
                               weight_field=None, units=units, extrema=extrema)

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

.. _eps-writer:

Publication-ready Figures
-------------------------

While the routines above give a convienent method to inspect and
visualize your data, publishers often require figures to be in PDF or
EPS format.  While the matplotlib supports vector graphics and image
compression in PDF formats, it does not support compression in EPS
formats.  The :class:`~yt.visualization.eps_writer.DualEPS` module
provides an interface with the `PyX <http://pyx.sourceforge.net/>`_,
which is a Python abstraction of the PostScript drawing model with a
LaTeX interface.  It is optimal for publications to provide figures
with vector graphics to avoid rasterization of the lines and text,
along with compression to produce figures that do not have a large
filesize.

.. note::
   PyX must be installed, which can be accomplished either manually
   with ``pip install pyx`` or with the install script by setting
   ``INST_PYX=1``.

This module can take any of the plots mentioned above and create an
EPS or PDF figure.  For example,

.. code-block:: python

    import yt.visualization.eps_writer as eps
    slc = yt.SlicePlot(ds, 'z', 'density')
    slc.set_width(25, 'kpc')
    eps_fig = eps.single_plot(slc)
    eps_fig.save_fig('zoom', format='eps')
    eps_fig.save_fig('zoom-pdf', format='pdf')

The ``eps_fig`` object exposes all of the low-level functionality of
``PyX`` for further customization (see the `PyX documentation
<http://pyx.sourceforge.net/manual/index.html>`_).  There are a few
convenience routines in ``eps_writer``, such as drawing a circle,

.. code-block:: python

    eps_fig.circle(radius=0.2, loc=(0.5,0.5))
    eps_fig.sav_fig('zoom-circle', format='eps')

with a radius of 0.2 at a center of (0.5, 0.5), both of which are in
units of the figure's field of view.  The
:func:`~yt.visualization.eps_writer.multiplot_yt` routine also
provides a convenient method to produce multi-panel figures
from a PlotWindow.  For example,

.. code-block:: python

    import yt
    import yt.visualization.eps_writer as eps
   
    slc = yt.SlicePlot(ds, 'z', ['density', 'temperature', 'pressure',
                       'velocity_magnitude'])
    slc.set_width(25, 'kpc')
    eps_fig = eps.multiplot_yt(2, 2, slc, bare_axes=True)
    eps_fig.scale_line(0.2, '5 kpc')
    eps_fig.save_fig('multi', format='eps')

will produce a 2x2 panel figure with a scale bar indicating 5 kpc.
The routine will try its best to place the colorbars in the optimal
margin, but it can be overridden by providing the keyword
``cb_location`` with a dict of either ``right, left, top, bottom``
with the fields as the keys.

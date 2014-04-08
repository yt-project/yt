
.. _how-to-make-plots:

How to Make Plots
=================

In this section we explain how to use ``yt`` to create visualizations
of simulation data, derived fields, and the data produced by ``yt``
analysis objects.  For details about the data extraction and
algorithms used to produce the image and analysis data, please see the
``yt`` `method paper
<http://adsabs.harvard.edu/abs/2011ApJS..192....9T>`_.  There are also
many example scripts in :ref:`cookbook`.

The :class:`~yt.visualization.plot_window.PlotWindow` interface is useful for
taking a quick look at simulation outputs.  Simple mechanisms exist for making 
plots of slices, projections, 1D profiles, and 2D profiles (phase plots), all of 
which are described below.

.. _simple-inspection:

Visual Inspection
-----------------

If you need to take a quick look at a single simulation output, ``yt``
provides the ``PlotWindow`` interface for generating annotated 2D
visualizations of simulation data.  You can create a ``PlotWindow`` plot by
supplying a parameter file, a list of fields to plot, and a plot center to
create a :class:`~yt.visualization.plot_window.SlicePlot`,
:class:`~yt.visualization.plot_window.ProjectionPlot`, or
:class:`~yt.visualization.plot_window.OffAxisSlicePlot`.

Plot objects use ``yt`` data objects to extract the maximum resolution
data available to render a 2D image of a field. Whenever a
two-dimensional image is created, the plotting object first obtains
the necessary data at the *highest resolution*.  Every time an image
is requested of it -- for instance, when the width or field is changed
-- this high-resolution data is then pixelized and placed in a buffer
of fixed size. This is accomplished behind the scenes using
:class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`
``PlotWindow`` expose the underlying matplotlib ``figure`` and
``axes`` objects, making it easy to customize your plots and 
add new annotations.

.. _slice-plots:

Slice Plots
-----------

The quickest way to plot a slice of a field through your data is to use
:class:`~yt.visualization.plot_window.SlicePlot`.  Say we want to visualize a
slice through the Density field along the z-axis centered on the center of the
simulation box in a simulation dataset we've opened and stored in the parameter
file object ``pf``.  This can be accomplished with the following command:

.. code-block:: python

   >>> slc = SlicePlot(pf, 'z', 'density')
   >>> slc.save()

These two commands will create a slice object and store it in a variable we've
called ``slc``.  We then call the ``save()`` function that is associated with
the slice object.  This automatically saves the plot in png image format with an
automatically generated filename.  If you don't want the slice object to stick
around, you can accomplish the same thing in one line:

.. code-block:: python
   
   >>> SlicePlot(pf, 'z', 'density').save()

It's nice to keep the slice object around if you want to modify the plot.  By
default, the plot width will be set to the size of the simulation box.  To zoom
in by a factor of ten, you can call the zoom function attached to the slice
object:

.. code-block:: python

   >>> slc = SlicePlot(pf, 'z', 'density')
   >>> slc.zoom(10)
   >>> slc.save('zoom')

This will save a new plot to disk with a different filename - prepended with
'zoom' instead of the name of the parameter file. If you want to set the width
manually, you can do that as well. For example, the following sequence of
commands will create a slice, set the width of the plot to 10 kiloparsecs, and
save it to disk.

.. code-block:: python

   >>> slc = SlicePlot(pf, 'z', 'density')
   >>> slc.set_width((10,'kpc'))
   >>> slc.save('10kpc')

The SlicePlot also optionally accepts the coordinate to center the plot on and
the width of the plot:

.. code-block:: python

   >>> SlicePlot(pf, 'z', 'density', center=[0.2, 0.3, 0.8], 
   ...           width = (10,'kpc')).save()

The center must be given in code units.  Optionally, you can supply 'c' or 'm'
for the center.  These two choices will center the plot on the center of the
simulation box and the coordinate of the maximum density cell, respectively.

Here is an example that combines all of the options we just discussed.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', center=[0.5, 0.5, 0.5], width=(20,'kpc'))
   slc.save()

The above example will display an annotated plot of a slice of the
Density field in a 20 kpc square window centered on the coordinate
(0.5, 0.5, 0.5) in the x-y plane.  The axis to slice along is keyed to the
letter 'z', corresponding to the z-axis.  Finally, the image is saved to
a png file.

Conceptually, you can think of the SlicePlot as an adjustable window
into the data. For example:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'pressure', center='c')
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

A slice object can also add annotations like a title, an overlying
quiver plot, the location of grid boundaries, halo-finder annotations,
and many other annotations, including user-customizable annotations.
For example:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.annotate_grids()
   slc.save()

will plot the density field in a 10 kiloparsec slice through the
z-axis centered on the highest density point in the simulation domain.
Before saving the plot, the script annotates it with the grid
boundaries, which are drawn as thick black lines by default.

Annotations are described in :ref:`callbacks`.  See
:class:`~yt.visualization.plot_window.SlicePlot` for the full class
description.

.. _projection-plots:

Projection Plots
----------------

Using a fast adaptive projection, ``yt`` is able to quickly project
simulation data along the coordinate axes.

Projection plots are created by instantiating a
:class:`~yt.visualization.plot_window.ProjectionPlot` object.  For
example:

.. python-script::
 
   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   prj = ProjectionPlot(pf, 2, 'density', width=(25, 'kpc'), 
                        weight_field=None)
   prj.save()

will create a projection of Density field along the x axis, plot it,
and then save it to a png image file.  The projection is only carried
out to level 2 of the AMR index and no weighting is applied.

Like :ref:`slice-plots`, annotations and modifications can be applied
after creating the ``ProjectionPlot`` object.  Annotations are
described in :ref:`callbacks`.  See
:class:`~yt.visualization.plot_window.ProjectionPlot` for the full
class description.

.. _off-axis-slices:

Off Axis Slice Plots
--------------------

Off axis slice plots can be generated in much the same way as
grid-aligned slices.  Off axis slices use
:class:`~yt.data_objects.data_containers.AMRCuttingPlaneBase` to slice
through simulation domains at an arbitrary oblique angle.  A
:class:`~yt.visualization.plot_window.OffAxisSlicePlot` can be
instantiated by specifying a parameter file, the normal to the cutting
plane, and the name of the fields to plot.  For example:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1,1,0] # vector normal to cutting plane
   north_vector = [-1,1,0]
   cut = OffAxisSlicePlot(pf, L, 'density', width=(25, 'kpc'), 
                          north_vector=north_vector)
   cut.save()

creates an off-axis slice in the plane perpendicular to ``L``,
oriented such that ``north_vector`` is the up direction.  If ``L`` and
``north_vector`` are not perpendicular, the component of
``north_vector`` perpendicular to ``L`` is used. Like
:ref:`slice-plots`, annotations and modifications can be applied after
creating the ``OffAxisSlicePlot`` object.  Annotations are described
in :ref:`callbacks`.  See
:class:`~yt.visualization.plot_window.OffAxisSlicePlot` for the full
class description.

.. _off-axis-projections:

Off Axis Projection Plots
-------------------------

Off axis projection plots .  Internally, off axis projections are
created using :ref:`the-camera-interface` by applying the
:class:`~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction`.
In this use case, the volume renderer casts a set of plane
parallel rays, one for each pixel in the image.  The data values
along each ray are summed, creating the final image buffer.

.. _off-axis-projection-function:

To avoid manually creating a camera and setting the transfer
function, yt provides the :func:`~yt.visualization.volume_rendering.camera.off-axis-projection`
function, which wraps the camera interface to create an off axis
projection image buffer.  These images can be saved to disk or
used in custom plots.  This snippet creates an off axis
projection through a simulation.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1,1,0] # vector normal to cutting plane
   north_vector = [-1,1,0]
   W = [0.02, 0.02, 0.02]
   c = [0.5, 0.5, 0.5]
   N = 512
   image = off_axis_projection(pf, c, L, W, N, "density")
   write_image(np.log10(image), "%s_offaxis_projection.png" % pf)

Here, ``W`` is the width of the projection in the x, y, *and* z
directions.

One can also generate generate annotated off axis projections
using
:class:`~yt.visualization.plot_window.OffAxisProjectionPlot`. These
plots can be created in much the same way as an
``OffAxisSlicePlot``, requiring only an open dataset, a direction
to project along, and a field to project.  For example:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   L = [1,1,0] # vector normal to cutting plane
   north_vector = [-1,1,0]
   prj = OffAxisProjectionPlot(pf,L,'density',width=(25, 'kpc'), 
                               north_vector=north_vector)
   prj.save()

OffAxisProjectionPlots can also be created with a number of
keyword arguments, as described in the `api reference`__ for the
class initializer.

__ :class:`~yt.visualization.plot_window.OffAxisProjectionPlot`

Plot Customization
------------------

You can customize each of the four plot types above in identical ways.  We'll go
over each of the customizations methos below.  For each of the examples below we
will modify the following plot.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.save()

Panning and zooming
~~~~~~~~~~~~~~~~~~~

There are three methods to dynamically pan around the data.  

:class:`~yt.visualization.plot_window.SlicePlot.pan` accepts x and y deltas in code
units.

.. python-script::

   from yt.mods import *
   from yt.units import kpc
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.pan((2*kpc, 2*kpc))
   slc.save()

:class:`~yt.visualization.plot_window.SlicePlot.pan_rel` accepts deltas in units relative
to the field of view of the plot.  

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.pan_rel((0.1, -0.1))
   slc.save()

:class:`~yt.visualization.plot_window.SlicePlot.zoom` accepts a factor to zoom in by.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.zoom(2)
   slc.save()

Set axes units
~~~~~~~~~~~~~~

:class:`~yt.visualization.plot_window.SlicePlot.set_axes_unit` allows the customization of
the axes unit labels.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_axes_unit('Mpc')
   slc.save()

Set the plot center
~~~~~~~~~~~~~~~~~~~

The :class:`~yt.visualization.plot_window.SlicePlot.set_center` function accepts a new
center for the plot, in code units.  New centers must be two element tuples.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_center((0.5, 0.5))
   slc.save()

Fonts
~~~~~

:class:`~yt.visualization.plot_window.SlicePlot.set_font` allows font costomization.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_font({'family': 'sans-serif', 'style': 'italic','weight': 'bold', 'size': 24})
   slc.save()

Colormaps
~~~~~~~~~

Each of these functions accept two arguments.  In all cases the first argument
is a field name.  This makes it possible to use different custom colormaps for
different fields tracked by the plot object.

To change the colormap for the plot, call the
:class:`~yt.visualization.plot_window.SlicePlot.set_cmap` function.  Use any of the
colormaps listed in the :ref:`colormaps` section.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_cmap('density', 'RdBu_r')
   slc.save()

The :class:`~yt.visualization.plot_window.SlicePlot.set_log` function accepts a field name
and a boolean.  If the boolean is :code:`True`, the colormap for the field will
be log scaled.  If it is `False` the colormap will be linear.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_log('density', False)
   slc.save()

Lastly, the :class:`~yt.visualization.plot_window.SlicePlot.set_zlim` function makes it
possible to set a custom colormap range.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_zlim('density', 1e-30, 1e-25)
   slc.save()

Set the size of the plot
~~~~~~~~~~~~~~~~~~~~~~~~

To set the size of the plot, use the
:class:`~yt.visualization.plot_window.SlicePlot.set_window_size` function.  The argument
is the size of the longest edge of the plot in inches.  View the full resolution
image to see the difference more clearly.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_window_size(10)
   slc.save()

To change the resolution of the image, call the
:class:`~yt.visualization.plot_window.SlicePlot.set_buff_size` function.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'))
   slc.set_buff_size(1600)
   slc.save()

.. _how-to-make-1d-profiles:

1D Profile Plots
----------------

1D profiles are used to calculated the average or the sum of a given quantity
with respect to a second quantity.  This means the "average density as a
function of radius" or "the total mass within a given set of density bins."
When created, they default to the average: in fact, they default to the average
as weighted by the total cell mass.  However, this can be modified to take
either the total value or the average with respect to a different quantity.

Profiles operate on data objects; they will take the entire data contained in a
sphere, a prism, an extracted region and so on, and they will calculate and use
that as input to their calculation.  To make a 1D profile plot, create a
(:class:`~yt.visualization.profile_plotter.ProfilePlot`) object, supplying the 
data object, the field for binning, and a list of fields to be profiled.

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_galaxy = pf.disk([0.5, 0.5, 0.5], [0.0, 0.0, 1.0], 0.01, 0.003)
   plot = ProfilePlot(my_galaxy, "density", ["temperature"])
   plot.save()

This will create a :class:`yt.data_objects.data_containers.AMRCylinderBase`
centered at [0.3, 0.5, 0.8], with a normal vector of [0.4, 0.5, 0.1], radius of
0.01 and height of 0.001 and will then make a plot of the average (as a 
function of the cell mass) temperature as a function of density.

As another example, we create a sphere of radius 100 pc and plot total mass 
in every equally-spaced temperature bin/

We could also have allowed the plot collection to create a sphere for us, as
well.  For instance:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = pf.sphere([0.5, 0.5, 0.5], (100, "kpc"))
   plot = ProfilePlot(my_sphere, "temperature", ["cell_mass"],
                      weight_field=None)
   plot.save()

Note that because we have specified the weighting field to be none, it operates 
as a local-bin accumulator.  We can also accumulate along the x-axis by setting 
the **accumulation** keyword argument to True, which is useful for plots of 
enclosed mass.

You can also access the data generated by profiles directly, which can be
useful for overplotting average quantities on top of phase plots, or for
exporting and plotting multiple profiles simultaneously from a time series.
The ``profiles`` attribute contains a list of all profiles that have been 
made.  For each item in the list, the x field data can be accessed with ``x``.  
The profiled fields can be accessed from the dictionary ``field_data``.

.. code-block:: python

   plot = ProfilePlot(my_sphere, "temperature", ["cell_mass"],
                      weight_field=None)
   # print the x field
   print plot.profiles[-1].x
   # print the profiled temperature field
   print plot.profiles[-1].field_data["temperature"]

Other options, such as the number of bins, are also configurable. See the 
documentation for 
The number of bins and other options and tweaks are 
available for these methods.  See the documentation for 
:class:`~yt.visualization.profile_plotter.ProfilePlot`
for more information.

Overplotting Multiple 1D Profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is often desirable to overplot multiple 1D profile to show evolution 
with time.  This is supported with the ``from_profiles`` class method.  
1D profiles are created with the :meth:`yt.data_objects.profiles.create_profile` 
method and then given to the ProfilePlot object.

.. python-script::

   from yt.mods import *

   # Create a time-series object.
   es = simulation("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo")
   es.get_time_series(redshifts=[5, 4, 3, 2, 1, 0])


   # Lists to hold profiles, labels, and plot specifications.
   profiles = []
   labels = []

   # Loop over each dataset in the time-series.
   for pf in es:
       # Create a data container to hold the whole dataset.
       ad = pf.h.all_data()
       # Create a 1d profile of density vs. temperature.
       profiles.append(create_profile(ad, ["temperature"], 
                                      fields=["cell_mass"],
                                      weight_field=None,
                                      accumulation=True))
       # Add labels
       labels.append("z = %.2f" % pf.current_redshift)

   # Create the profile plot from the list of profiles.
   plot = ProfilePlot.from_profiles(profiles, labels=labels)

   # Save the image.
   plot.save()

Altering Line Properties
~~~~~~~~~~~~~~~~~~~~~~~~

Line properties for any and all of the profiles can be changed with the 
``set_line_property`` function.  The two arguments given are the line 
property and desired value.

.. code-block:: python

   >>> plot.set_line_property("linestyle", "--")

With no additional arguments, all of the lines plotted will be altered.  To 
change the property of a single line, give also the index of the profile.

.. code-block:: python

   >>> # change only the first line
   >>> plot.set_line_property("linestyle", "--", 0)

.. _how-to-make-2d-profiles:

2D Phase Plots
--------------

2D phase plots function in much the same was as 1D phase plots, but with a 
:class:`~yt.visualization.profile_plotter.PhasePlot` object.  Much like 1D profiles, 
2D profiles (phase plots) are best thought of as plotting a distribution of points, 
either taking the average or the accumulation in a bin.  For example, to generate a 
2D distribution of mass enclosed in density and temperature bins, you can do:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = pf.sphere("c", (50, "kpc"))
   plot = PhasePlot(my_sphere, "density", "temperature", ["cell_mass"],
                    weight_field=None)
   plot.save()

If you would rather see the average value of a field as a function of two other
fields, you can neglect supplying the *weight* parameter.  This would look
something like:

.. python-script::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   my_sphere = pf.sphere("c", (50, "kpc"))
   plot = PhasePlot(my_sphere, "density", "temperature", ["H_fraction"],
                    weight_field="cell_mass")
   plot.save()

Probability Distribution Functions and Accumulation
---------------------------------------------------

Both 1D and 2D profiles which show the total of amount of some field, such as mass, 
in a bin (done by setting the **weight_field** keyword to None) can be turned into 
probability distribution functions (PDFs) by setting the **fractional** keyword to 
True.  When set to True, the value in each bin is divided by the sum total from all 
bins.  These can be turned into cumulative distribution functions (CDFs) by setting 
the **accumulation** keyword to True.  This will make is so that the value in any bin 
N is the cumulative sum of all bins from 0 to N.  The direction of the summation can be 
rversed by setting **accumulation** to -True.  For PhasePlots, the accumulation can be 
set independently for each axis by setting **accumulation** to a list of True/-True/False 
values.

.. _interactive-plotting:

Interactive Plotting
--------------------

The best way to interactively plot data is through the IPython notebook.  Many
detailed tutorials on using the IPython notebook can be found at
http://ipython.org/presentation.html , but the simplest way to use it is to
type:

.. code-block:: bash

   yt notebook

at the command line.  This will prompt you for a password (so that if you're on
a shared user machine no one else can pretend to be you!) and then spawn an
IPython notebook you can connect to.

If you want to see yt plots inline inside your notebook, you need only create a
plot and then call ``.show()``:

.. notebook-cell::

   from yt.mods import *
   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   p = ProjectionPlot(pf, "x", "density", center='m', width=(10,'kpc'),
                      weight_field='density')
   p.show()

The image will appear inline.

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

   >>> import yt.visualization.eps_writer as eps
   >>> slc = SlicePlot(pf, 'z', 'density')
   >>> slc.set_width(25, 'kpc')
   >>> eps_fig = eps.single_plot(slc)
   >>> eps_fig.save_fig('zoom', format='eps')
   >>> eps_fig.save_fig('zoom-pdf', format='pdf')

The ``eps_fig`` object exposes all of the low-level functionality of
``PyX`` for further customization (see the `PyX documentation
<http://pyx.sourceforge.net/manual/index.html>`_).  There are a few
convenience routines in ``eps_writer``, such as drawing a circle,

.. code-block:: python

   >>> eps_fig.circle(radius=0.2, loc=(0.5,0.5))
   >>> eps_fig.sav_fig('zoom-circle', format='eps')

with a radius of 0.2 at a center of (0.5, 0.5), both of which are in
units of the figure's field of view.  The
:class:`~yt.visualization.eps_writer.multiplot_yt` routine also
provides a convenient method to produce multi-panel figures
from a PlotWindow.  For example,

.. code-block:: python

   >>> import yt.visualization.eps_writer as eps
   >>> slc = SlicePlot(pf, 'z', ['density', 'temperature', 'Pressure',
                       'VelocityMagnitude'])
   >>> slc.set_width(25, 'kpc')
   >>> eps_fig = eps.multiplot_yt(2, 2, slc, bare_axes=True)
   >>> eps_fig.scale_line(0.2, '5 kpc')
   >>> eps_fig.save_fig('multi', format='eps')

will produce a 2x2 panel figure with a scale bar indicating 5 kpc.
The routine will try its best to place the colorbars in the optimal
margin, but it can be overridden by providing the keyword
``cb_location`` with a dict of either ``right, left, top, bottom``
with the fields as the keys.

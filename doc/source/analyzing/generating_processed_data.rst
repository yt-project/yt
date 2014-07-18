.. _generating-processed-data:

Generating Processed Data
=========================

Although ``yt`` provides a number of built-in visualization methods that can
process data and construct from that plots, it is often useful to generate the
data by hand and construct plots which can then be combined with other plots,
modified in some way, or even (gasp) created and modified in some other tool or
program.

.. _generating-2d-image-arrays:

2D Image Arrays
---------------

When making a slice, a projection or an oblique slice in yt, the resultant
:class:`~yt.data_objects.data_containers.AMR2DData` object is created and
contains flattened arrays of the finest available data.  This means a set of
arrays for the x, y, (possibly z), dx, dy, (possibly dz) and data values, for
every point that constitutes the object.


This presents something of a challenge for visualization, as it will require
the transformation of a variable mesh of points consisting of positions and
sizes into a fixed-size array that appears like an image.  This process is that
of pixelization, which ``yt`` handles transparently internally.  You can access
this functionality by constructing a
:class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer` (or 
:class:`~yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer`) and
supplying to it your :class:`~yt.data_objects.data_containers.AMR2DData`
object, as well as some information about how you want the final image to look.
You can specify both the bounds of the image (in the appropriate x-y plane) and
the resolution of the output image.  You can then have ``yt`` pixelize any
field you like.

To create :class:`~yt.data_objects.data_containers.AMR2DData` objects, you can
access them as described in :ref:`using-objects`, specifically the section
:ref:`available-objects`.  Here is an example of how to window into a slice 
of resolution(512, 512) with bounds of (0.3, 0.5) and (0.6, 0.8).  The next
step is to generate the actual 2D image array, which is accomplished by
accessing the desired field.

.. code-block:: python

   sl = ds.slice(0, 0.5)
   frb = FixedResolutionBuffer(sl, (0.3, 0.5, 0.6, 0.8), (512, 512))
   my_image = frb["density"]

This resultant array can be saved out to disk or visualized using a
hand-constructed Matplotlib image, for instance using
:func:`~matplotlib.pyplot.imshow`.

.. _generating-profiles-and-histograms:

Profiles and Histograms
-----------------------

Profiles and histograms can also be generated using the
:class:`~yt.visualization.profile_plotter.ProfilePlot` and 
:class:`~yt.visualization.profile_plotter.PhasePlot` functions 
(described in :ref:`how-to-make-1d-profiles` and
:ref:`how-to-make-2d-profiles`).  These generate profiles transparently, but the
objects they handle and create can be handled manually, as well, for more
control and access.  For instance, if you wanted to plot a time series of the
evolution of a profile, or if you wanted to handle the fields in some way to
calculate an accretion rate or some modified version of the resultant
histogram.  For full documentation, see the API reference for
:class:`~yt.data_objects.profiles.BinnedProfile1D`,
:class:`~yt.data_objects.profiles.BinnedProfile2D`, or
:class:`~yt.data_objects.profiles.BinnedProfile3D`.

Profile objects can be created from any data object (see :ref:`using-objects`,
specifically the section :ref:`available-objects` for more information) and are
best thought of as distribution calculations.  They can either sum up or
average one quantity with respect to one or more other quantities, and they do
this over all the data contained in their source object.

To generate a profile, you need to supply the limits of the distribution for
each variable along which you are distributing (i.e., the x- and y-axes for 2D
profiles, but only the x-axis for 1D profiles) as well as the number of bins
into which you want the values distributed.  Often these are the least
straightforward pieces of information; the usage of derived quantities,
specifically ``Extrema``, can help with this.  (See :ref:`derived-quantities`
for more information on this.)  Once you have created the profile object, you
can add fields to it either one at a time or multiple simultaneously.  If you
supply a weighting field, the average will be taken.  Otherwise, if the weight
field is set to ``None``, only an accumulation inside a bin will be performed.
Note that by default the weight field is ``CellMassMsun``!

For instance, to create a sphere at (0.3, 0.6, 0.4) and then take the 1D
average distribution of fields with respect to Density, you would first
construct your profile.  Then you would add fields to it; for instance, we can
add ``CellMassMsun`` in an unweighted fashion to get the total mass in each
bin.  Then we add ``Temperature`` with the default weighting to get the
average value in each bin.  Here's an example, where we have used our knowledge
of the bounds of density in advance to set up the profile.

.. code-block:: python

   source = ds.sphere( (0.3, 0.6, 0.4), 1.0/ds['pc'])
   profile = BinnedProfile1D(source, 128, "density", 1e-24, 1e-10)
   profile.add_fields("cell_mass", weight = None)
   profile.add_fields("temperature")

At this point, we can access the fields ``CellMassMsun`` and ``Temperature``
from the ``profile`` object, which are returned as 1D arrays.

.. code-block:: python

   print profile["cell_mass"]
   print profile["temperature"]

The field ``UsedBins`` is also included, which is ``True`` wherever values have
been added.  This is primarily used for 2D profiles, where many of the bins may
be empty and need to be masked.  Note also that the bins used to generate the
profiles, in this case ``Density``, are also defined to allow for x-y plots.

One of the more interesting techniques that is enabled with this approach is
the generation of 1D profiles that correspond to 2D profiles.  For instance, a
phase plot that shows the distribution of mass in the density-temperature
plane, with the average temperature overplotted.

To generate a 2D profile, the interface is broadly the same except with a few
additional parameters for the second field along which values will be
distributed.  Here we are also distributing values along temperature, and then
calculating the mass in each (2D) bin.

.. code-block:: python

   source = ds.sphere( (0.3, 0.6, 0.4), 1.0/ds['pc'])
   prof2d = BinnedProfile2D(source, 128, "density", 1e-24, 1e-10, True,
                                    128, "temperature", 10, 10000, True)
   prof2d.add_fields("cell_mass", weight = None)

Note that at this point we can use :func:`~matplotlib.pyplot.pcolormesh` to
plot the ``prof2d["cell_mass"]`` value, and even overplot the value of
``profile["temperature"]`` to show the average value in every density bin.
Note that you will likely have to mask out the zero values using the
``prof2d["UsedBins"]`` field.  Profiles can also be calculated in
three-dimensions, with a similar extension of the calling function.

.. _generating-line-queries:

Calculating the Variance of Profiled Fields
+++++++++++++++++++++++++++++++++++++++++++

See :ref:`cookbook-profile-variance` for an example of the following.  
When calculating average 1D and 2D profiles (when the *weight* keyword is not 
None), the variance within each bin is calculated automatically.  A practical 
application for this would be calculating velocity dispersion by profiling the 
average velocity magnitude.  The variance values for 1D and 2D profiles are 
accessible as the name of the profiled field followed by ``_std``.  For the 
above examples, this is done with

.. code-block:: python

   print profile["Temperature_std"]
   print prof2d["Temperature_std"]

Line Queries and Planar Integrals
---------------------------------

To calculate the values along a line connecting two points in a simulation, you
can use the object :class:`~yt.data_objects.data_containers.AMRRayBase`,
accessible as the ``ray`` property on a index.  (See :ref:`using-objects`
for more information on this.)  To do so, you can supply two points and access
fields within the returned object.  For instance, this code will generate a ray
between the points (0.3, 0.5, 0.9) and (0.1, 0.8, 0.5) and examine the density
along that ray:

.. code-block:: python

   ray = ds.ray(  (0.3, 0.5, 0.9), (0.1, 0.8, 0.5) )
   print ray["density"]

The points are ordered, but the ray is also traversing cells of varying length,
as well as taking a varying distance to cross each cell.  To determine the
distance traveled by the ray within each cell (for instance, for integration)
the field ``dt`` is available; this field will sum to 1.0, as the ray's path
will be normalized to 1.0, independent of how far it travels through the domain.
To determine the value of ``t`` at which the ray enters each cell, the field
``t`` is available.  For instance:

.. code-block:: python

   print ray['dts'].sum()
   print ray['t']

These can be used as inputs to, for instance, the Matplotlib function
:func:`~matplotlib.pyplot.plot`, or they can be saved to disk.

The volume rendering functionality in yt can also be used to calculate
off-axis plane integrals, using the
:class:`~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction`
in a manner similar to that described in :ref:`volume_rendering`.

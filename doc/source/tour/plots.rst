Plots
=====

Through the plotting interface, you can have ``yt`` automatically generate many
of the analysis objects available to you!

The primary plotting interface is through a :class:`PlotCollection`
instantiated with a given parameter file and (optionally) a center.  See
:ref:`quick_making_plots` for a brief example of how to generate a
:class:`PlotCollection`.

Two-Dimensional Images
----------------------

Whenever a two-dimensional image is created, the plotting object first
obtains the necessary data at the *highest resolution*.  Every time an image is
requested of it -- for instance, when the width or field is changed -- this
high-resolution data is then pixelized and placed in a buffer of fixed size.

Slices are axially-aligned images of data selected at a fixed point on an axis;
these are the fastest type of two-dimensional image, as only the correct
coordinate data is read from disk and then plotted.

Cutting planes are oblique slices, aligned with a given normal vector.  These
can be used for face-on images of disks and other objects, as well as a
rotational slices.  They work just like slices in other ways, but they tend to
be a bit slower.

Projections are closer in style to profiles than slices.  They can exist either
as a summation of the data along every possible ray through the simulation, or
an average value along every possible ray.  If a *weight_field* is provided,
then the data returned is an average; typically you will want to weight with
``Density``.  If you do not supply a *weight_field* then the returned data is a
column sum.  These fields are stored in between invocations -- this allows for
speedier access to a relatively slow process!

.. _profiles_and_phase_plots:

Profiles and Phase Plots
------------------------

Profiles and phase plots provide identical API to the generation of profiles
themselves, but with a couple convenience interfaces.  You can have the plot
collection generate a sphere automatically for either one:

.. code-block:: python

   pc.add_phase_sphere(100.0, 'au', ["Density", "Temperature", "CellMassMsun"],
                       weight = None)

This will generate a sphere, a phase plot, and then return to you the plot
object.

Interactive Plotting
--------------------

Thanks to the pylab interface in Matplotlib, we have an interactive plot
collection available for usage within ``IPython``.  Instead of
:class:`PlotCollection`, use :class:`PlotCollectionInteractive` -- this will
generate automatically updating GUI windows with the plots inside them.

Callbacks
---------

Callbacks are means of adding things on top of existing plots -- like vectors,
overplotted lines, and so on and so forth.  They have to be added to the plot
objects themselves, rather than the :class:`PlotCollection`.  You can add them like so:

.. code-block:: python

   p = pc.add_slice("Density", 0)
   p.add_callback(GridBoundaryCallback())

Each Callback has to be instantiated, and then added.  You can also access the
plot objects inside the PlotCollection directly:

.. code-block:: python

   pc.add_slice("Density", 0)
   pc.plots[-1].add_callback(GridBoundaryCallback())

Note that if you are plotting interactively, the PlotCollection will need to
have ``redraw`` called on it.

For more information about Callbacks, see the :ref:`API reference <callbacks>` .

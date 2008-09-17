.. index::
   single: GUI
   module: reason
   module: yt.reason

Using Reason
============

Reason is the GUI for yt.  It's developed in wxPython, and all attempts have
been made to expose a susbtantial API to enable interactive plot addition, as
well as a command-line for examining raw data.

Starting Up Reason
------------------

If you've downloaded the binary package for OSX, you can simply double-click
the icon.  Otherwise, if you have `wxPython <http://www.wxpython.org/>`_
installed, and yt is installed globally, you should be able to simply type
``reason`` to initiate it.

When reason starts up, it will execute all the code in your plugin file, so any
fields you've defined there will be defined in the dialog boxes.

Go ahead and open up a dataset now, using the file menu's "Open Hierarchy"
option.  It'll get added to the tree on the left.  All of its subsequent
derived data objects -- slices, projections, spheres, etc -- will be added
there.

.. index:: projection, slice

Making Slices and Projections
-----------------------------

By right clicking on an output file, you can generate a slice.  The software
will automatically add a slice for all three axes; with projections, you're
given more control over the field, weight and axes.

Each created plot will reside in its own tabs; these tabs are linked in their
display of fields, their width, and their center.  They can be rearranged by
clicking and dragging -- for instance, multiple panes can be viewed at a single
time by clicking and dragging to the region you'd like to split into.

These plots can be recentered by getting the context menu on the plot image and
choosing to either center on the location of the click or to center on the
location of maximum density.

.. note:: On some platforms, the right-click on a plot in wxPython is
   non-functional; on these platforms, double-click to pop up the context menu.

Give it a shot; re-center somewhere else.  All of the plots generated at a
given time will be linked, so they should all re-center on wherever you told
them to.

.. note:: Even though a projection loses some information, by re-centering on
   two of the three axes, you can uniquely determine a new center.

.. index:: sphere, EnzoSphere

Obtaining Spheres
-----------------

So now we've got a slice through our data, and we've dealt with recentering,
changing the width and so on.  Now let's try something more advanced.

You can generate a sphere by two methods.  The first is to bring up the context
menu on a plot (either right- or double-click it) and then requesting one.
This will be centered at the current "plot center" and you can specify the
width exactly.

Alternatively, you can hold down the *shift* key and click once on the center
and once on the outer edge.  This will generate a sphere using those two pieces
of information as the defining characteristics.

This sphere will appear in the object tree as a child of the current plot.

.. index:: phase plot

Making Phase and Profile Plots
------------------------------

Now that you have a sphere, you are able to make phase and profile plots from
it.  Right click on it in the object tree and choose to make either a phase
(two-dimensional) or profile (one-dimensional) plot.

Each of these choices will present you with a constructor window, where you can
enter the appropriate parameters to define the plot.  Each of these options is
reflected in the API for the :class:`BinnedProfile1D` and :class:`BinnedProfile2D`
classes respectively.

.. index:: cutting plane

Making Cutting Planes
---------------------

Within the GUI, once you have selected a sphere, you can generate a cutting
plane based on the mass-weighted angular momentum vector of that sphere, and
centered at its center.  You can access this by right-clicking on a sphere and
choosing "Cutting Plane."

Additionally, you can overlay a vector plot of the velocity vectors on this
plot via the context menu.

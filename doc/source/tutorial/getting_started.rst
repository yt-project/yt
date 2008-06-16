Getting Started
===============

.. index:: Starting yt

Getting into yt
---------------

For the purposes of this tutorial, I'm going to assume that python is
installed, the prerequisites are installed and working, and that you know how
to launch a python interpreter.  So first off, let's start it up.

.. code-block:: bash

   $ python2.5

We start out by getting some stuff into our local namespace.
The import command is how an external module gets loaded.  By default,
it shows up as the module name itself.  But you can also import it as
something else, which makes things easier...

Lagos is our data-handler.  If you need to touch data, lagos will be your
friend. ::

   >>> import yt.lagos as lagos

The :mod:`yt.lagos` module is now in our namespace and called :mod:`lagos`; you
can import it as whatever you like, but we'll stick with that for the remainder
of the tutorial.  Let's move on to opening an output.

.. index:: hierarchy

The Hierarchy
-------------

We first instantiate a StaticOutput, which is a single disk-written timestep of
the simulation  This is a pretty simple process - :mod:`lagos`
will grab the parameter file that you give it, parse it, set up some handy
units, and then return control to you.::

   >>> a = lagos.EnzoStaticOutput("my_data0001.dir/my_data0001")


But what good is just a parameter file?  What we really want is to be able to
manipulate the grids, and all the fields stored in those grids.  Fortunately,
we can access the hierarchy as either 'hierarchy' or 'h'.

The 'h' property will instantiate an :class:`EnzoHierarchy` object, which
posseses information about the grids in a given output, as well as their
spatial relationships to each other.  Furthermore, :mod:`yt` will store
commonly-used and expensive-to-generate data in an HDF5 file associated with
that object -- things like projections, 3D phase profiles, and so on.

All of the data-handling we want to do will be mediated by the hierarchy file;
so let's start out by finding the location of the maximum Density.
Note that here we use an inline-unpacking -- find_max returns two values,
packed into a tuple, so we manually unpack them right here.::

   >>> v,c = a.h.find_max("Density")
   >>> print v
   >>> print c

.. index:: Grid, 3ddata;grid

Now we know where the maximum density is located, and what it is.  This is a
pretty simple process, and we will use it later to center spheres and boxes and
disks on the point of maximum density.


The Grids
---------

The hierarchy object has some interesting methods that let us access the data
in different ways.  Let's start by looking at the grids.  ::

   >>> my_grid = a.h.grids[0] # The first grid in the hierarchy
   >>> print my_grid.dx, my_grid.ActiveDimensions # Some of the attributes


You might be asking yourself, how do I know what properties objects have?
Python gives us a convenient mechanism for determining this: ::

   >>> print dir(my_grid)

In general, a method or property that begins with '_' is not meant to be
called directly, and '__' means that it definitely shouldn't be called by an
external process.

Now let's take a look at accessing the data fields inherent to the grid
object.  ::

   >>> print my_grid["Density"]*my_grid["Temperature"]

What the code just did was access the grid files and then spit back the
results of the read operations.  So we can access the fields pretty easily
this way, but what about generating new fields from our existing fields?
Well, that's just as easy!  And yt comes with many derived fields
pre-defined, but from within a script you can also define your own.  In the
next section we'll take about the difference between a native field and a
derived field.

.. index:: derivedfield, indices

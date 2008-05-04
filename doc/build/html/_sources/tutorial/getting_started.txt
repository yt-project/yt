Getting Started
===============

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

The Hierarchy
-------------

We first instantiate a StaticOutput.  This is a pretty simple process - lagos
will grab the parameter file that you give it, parse it, set up some handy
units, and then return control to you.::

   >>> a = lagos.EnzoStaticOutput("/Users/matthewturk/Research/data/DataDir0036/DataDump0036")


But what good is just a parameter file?  What we really want is to be able to
manipulate the grids, and all the fields stored in those grids.  Fortunately,
we can access the hierarchy as either 'hierarchy' or 'h'.
Let's start out by finding the location of the maximum Density.
Note that here we use an inline-unpacking -- find_max returns two values,
packed into a tuple, so we manually unpack them right here.::

   >>> v,c = a.h.find_max("Density")
   >>> print v
   >>> print c

The Grids
---------

The hierarchy object has some interesting methods that let us access the data
in different ways.  Let's start by looking at the grids.  ::

   >>> my_grid = a.h.grids[0] # The root grid
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



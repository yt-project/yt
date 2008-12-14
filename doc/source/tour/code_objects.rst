Code Objects
============

Unfortunately, this name is a bit misleading, because in a sense, everything is
a code object!  But what it's meant to convey is that there are some objects
that are *native* to the AMR code that output them -- for Enzo, this would be
things like grid patches, the hierarchy, the parameter file, and a couple more
supplemental files that don't show up in analysis very often.

For the most part, ``yt`` will abstract these objects away -- there's no reason
the user (you & me!) should have to know anything about grids if they don't
*have* to.  But, at the same time, you should be able to dig as deep as you
want into the internals of the AMR code, and what it outputs.

Let's start with the structure of the data.  We can think about this in terms
of the *parameters* that describe the simulation, the *hierarchy* of grid
patches inside that simulation, and then the *grid patches* themselves.

Parameter Files
---------------

This is the most primitive and isolated means of examining data output from the
code, and it's also the jumping point for analyzing and plotting of any data
from a simulation -- this is your entry point inside ``yt``.  These particular
objects are meant to be very 'cheap' -- they're meant to be quick to
instantiate, so that you can, for instance, very rapidly identify which output
you want to examine in a script.  The class name is :class:`EnzoStaticOutput`.

For example, here's a snippet from my own research.  It uses a python module
called ``glob`` to do wildcard filename matching, and then it checks to see if
the time a given output was made at was greater than 200 million years, and if
so, it quits the loop, and then the script can process the file in whatever way.

.. code-block:: python

   from yt.mods import *
   import glob

   for fn in glob.glob("DataDump*"):
        pf = EnzoStaticOutput(fn)
        if pf["InitialTime"] * pf["years"] >= 2e8: break

This particular snippet also brings up the point that the parameter file
encompasses all of the units associated with a given output file.  It will do
its best to figure out how to convert between code quantities -- things like
time, energy, density, velocity -- and cgs quantities.

As a note, the parameter file actually has the ability to handle a couple other
code objects -- the cooling table and the rate table.  Right now the means of
processing this information is not terribly standardized, but please feel free
to experiment with the ``cool_rates`` and ``cool_rates`` objects.

Hierarchy
---------

With the hierarchy, we begin to delve a bit deeper into the workings of both
``yt`` and the AMR code.  This will parse the hiearchy file and instantiate all
of the grid objects inside python.  Unfortunately, with a large number of
grids, this can take some time.  (Efforts are under way to reduce this time,
but it's an unfortunate tradeoff between flexibility and speed.  Typically it
doesn't feel sluggish unless you are instantiating more than 50,000 grid
patches.)  This class is :class:`EnzoHierarchyType`, but you will only need to
access it as a property of an :class:`EnzoStaticOutput` called ``h``:

.. code-block:: python

   my_hierarchy = pf.h
   my_hierarchy.print_stats()

This will access the hierarchy, assign it to a local variable (not strictly
necessary!) and then print out statistics about the simulation -- the number of
levels, the number of grids on a level, the number of total cells, and some
information about the finest grid cell.  The hierarchy does have one very
important function called ``find_max`` -- this function will dig through the
top couple levels of a simulation and attempt to find the maximum point of a
given field, and then it returns the value and the center point to you.
(There's an example down below.)

Because the hierarchy contains all of the grids (in an array called ``grids``),
it can operate on the information about those grids very easily.  Additionally,
it acts as a means of generating the next level of objects, as well; in a
sense, the hierarchy acts as gatekeeper to the grid patches as well as to the
data in the simulation.

In the next section, we'll talk about some of these 'physical' objects, and how
you can use them to examine data in simulations.

Grids
-----

Inside a simulation, the data is first divided up into grid patches, which are
subdivided into cells, associated with each of which is a set of data points.
The grid patches may trace physical objects, but the distribution of grid
patches is governed not only by hydrodynamic structure but can also be
influenced by processor arrangement and topology.  As such, these are thought
of a code objects rather than physical objects.  However, despite that
categorization, the data returned to you when accessed through a grid is in cgs
units, as to allow better and more meaningful analysis.  This class is the
:class:`EnzoGridBase`.

You should not have to instantiate grid objects yourself -- they are
automatically created by the hierarchy, and stored in an array called
``grids``:

.. code-block:: python

   print my_hierarchy.grids[0].LeftEdge

Note that while the grids are indexed-by-one by Enzo in ``yt`` they are
indexed-by-zero in the grids array!

Usage of Code Objects
---------------------

For the most part, you should interact with the code objects in only a few
ways.  If ``yt`` has done its job correctly, you won't have to interact with
the grids individually or with the hierarchy except to create objects.
Typically you will do something like:

.. code-block:: python

   pf = EnzoStaticOutput("DataDump0019.dir/DataDump0019")
   v, c = pf.h.find_max("Density")
   sphere = pf.h.sphere(c, 100.0/pf['au'])

We now have a sphere (see the next section) centered at the point of maximum
density in the simulation.  At this point, we should be able to mostly
disregard the code objects, and operate exclusively on the sphere.

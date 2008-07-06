Advanced Objects
================

So we've seen simple objects already -- sphere, profiles, and regions.
There are two other types of objects that are fairly cool, as well, and
they behave the same way all the others do -- you can get any field from them,
profile them, all of that, right from the command line like you would a sphere
or a region.

.. index:: Cutting Plane, 2ddata;Cutting Plane

Cutting Plane
-------------

The first of these is the cutting plane.  Let's say you're doing some
collapse-based problem, where you end up with a vector that describes the
angular momentum.  By feeding this vector into the cutting plane object, you can
create an object that is sliced using that vector as the normal to the cutting
plane.  (For more information about how this is done, see the source code.)
An interface to this is available through raven.::

   >>> cutting_plane = a.h.cutting([0.1,0.3,0.4], [0.5,0.5,0.5])
   >>> pc.add_cutting_plane("Density", [0.1,0.3,0.4], center=[0.5,0.5,0.5])

These two statements create identical cutting-plane objects.  The latter
also prepares a plot to display the object, however.  Each of these objects
is uniquely determined by a normal vector and a point on the plane; we have
supplied the normal vector {0.1, 0.3, 0.4} and the center point of the system.
(Note that the normal vector need not be a unit vector, as it will be normalized
during the object instantiation.)  An appropriate 'up' vector is guessed at
and then used to define the transformed coordinate system.

.. index:: Extracted Regions, 3ddata;Extracted Regions

Extracted Regions
-----------------

One of the other interesting data types is mostly only useful for profiles in
one and two dimensions.  By whatever criteria you like, you can extract a subset
of any 3D data object.  What this means is that you can apply 'cuts' to the data --
whether it be by some tracer field or even just by a minimumt temperature --
and then receive a fully-function data object.

So for the next example, we will presuppose that we started out with two
galaxies, one of which was pre-populated with the BaryonField "GalaxyOne"
and the other with "GalaxyTwo".  As they collide, the gas will mix.  If we
wanted to profile all the gas that originated in one galaxy, we could create
a data object: ::

   >>> import numpy
   >>> gal1_region = a.h.sphere([0.25,0.25,0.25], 1.0/a["mpc"])
   >>> indices_we_want = numpy.where(gal1_region["GalaxyTwo"] > 1e-3)
   >>> gal1_mixed = gal1_region.extract_region(indices_we_want)

We now have a new data source that only contains the points of interest --
only the points where the contribution from the second galaxy is greater than
0.001, and we can then feed that as a data source into a Binned2DProfile.
Note that we used numpy here.  yt uses numpy behind the scenes for just about
everything, and we have been, too, when we've been using our data.  The 'where'
command creates a list of indices, which can then get fed into an array accessor
method, and that's the functionality we use here.

.. index:: GridCollection, 3ddata;GridCollection

Derived Quantities
==================

Associated with a given data object, one can imagine various quantities being
described.  Several of these are already included in yt, and more can be added.
To access a derived quantity, the property 'quantities' is available in a given
data object.  It acts as a dictionary, returning callable functions.::

   >>> L_vec = gal1_region.quantities["AngularMomentumVector"]()
   >>> com = gal1_region.quantities["CenterOfMass"]()
   >>> avg_T = gal1_region.quantities["WeightedAverageQuantity"]("Temperature", \
   ...                                                           "CellMassMsun")

All derived quantities accept the parameter :keyword:`lazy_reader`, which
tells the function to attempt to operate in a memory-conservative fashion --
reading and then flushing the data back to disk.  This is not possible for all
functions, notably the :func:`IsBound` function, which checks for
gravitationally-bound status.


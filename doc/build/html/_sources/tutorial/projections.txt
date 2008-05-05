Projections
===========

Projections are the line integrals of a given quantity along a given access.
These are useful in a wide variety of problems -- they can reveal morphologies
that slicing cannot, in the simplest case.

Simple Projections
------------------

There are several ways to generate projections in yt, and with the appropriate
software package installed (mpi4py) one can even generate them in parallel,
distributed across an arbitrary number of nodes.  However, I will deal here
with the simplest case -- generating a projection on a single processor.  Note that
because it's a time consuming task, yt will automatically store it in an HDF5
file so that you only have to do it once.

Much like cutting planes and slices, there are two ways to generate these --
one driven by the plotting engine raven, and the other as a direct call to
the hierarchy object.  ::

   >>> a.h.proj(0,"Temperature", weight_field="CellMass")
   >>> pc.add_projection("Temperature", 0, weight_field="CellMass")

The first call to project something with a given weight along a given axis
will cause that projection to be generated -- to full resolution -- and then
stored in an HDF5 file tied to the parameter file.  Then next time you ask
for it, it is already available.  Zooming and panning is accomplished very
simply.

Parallel Projections
--------------------



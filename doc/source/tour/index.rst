A Tour of yt
============

There's a lot inside yt.  This section is designed to give you an idea of
what's there, what you can do with it, and how to think about the yt
environment.

You will be served best by thinking about doing things with yt by thinking in
objects; we'll start by talking about how to access different types of objects
that are put on the disk by the AMR code, then we'll talk about how to select
subsections of a set of output based on physical quantities or other
characteristics, and then a bit about how to perform data reduction to
transform those objects into profiles, phase-space distributions, images, and
so on.  So this leads to a clear set of categories for different 'objects'
inside the yt-domain.

   Code Objects
     These are objects that are on the disk, that the AMR code knows about --
     things like grids, data dumps, the grid hierarchy and so on and so forth.
   Physical Objects
     These are objects like spheres, rectangular prisms, slices, and so on.
     These are collections of related data, and these associations are not
     necessarily associated with arrangement inside a code object.
   Reduced Objects
     These are objects where some set of data has been reduced into a more
     tractable format.  Histograms, 1-D profiles, averages and so on are all
     members of this category.
   Plots
     Plots stand somewhat outside this category, because the plotting interface
     accepts information about what you want to see, then goes out and fetches
     it as appropriate from the correct types of objects.

Contents:

.. toctree::
   :maxdepth: 2

   code_objects
   physical_objects
   reduced_objects
   plots

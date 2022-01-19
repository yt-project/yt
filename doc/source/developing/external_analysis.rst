.. _external-analysis-tools:

Using yt with External Analysis Tools
=====================================

yt can be used as a ``glue`` code between simulation data and other methods of
analyzing data.  Its facilities for understanding units, disk IO and data
selection set it up ideally to use other mechanisms for analyzing, processing
and visualizing data.

Calling External Python Codes
-----------------------------

Calling external Python codes very straightforward.  For instance, if you had a
Python code that accepted a set of structured meshes and then post-processed
them to apply radiative feedback, one could imagine calling it directly:

.. code-block:: python

   import radtrans

   import yt

   ds = yt.load("DD0010/DD0010")
   rt_grids = []

   for grid in ds.index.grids:
       rt_grid = radtrans.RegularBox(
           grid.LeftEdge,
           grid.RightEdge,
           grid["density"],
           grid["temperature"],
           grid["metallicity"],
       )
       rt_grids.append(rt_grid)
       grid.clear_data()

   radtrans.process(rt_grids)

Or if you wanted to run a population synthesis module on a set of star
particles (and you could fit them all into memory) it might look something like
this:

.. code-block:: python

   import pop_synthesis

   import yt

   ds = yt.load("DD0010/DD0010")
   ad = ds.all_data()
   star_masses = ad["StarMassMsun"]
   star_metals = ad["StarMetals"]

   pop_synthesis.CalculateSED(star_masses, star_metals)

If you have a code that's written in Python that you are having trouble getting
data into from yt, please feel encouraged to email the users list and we'll
help out.

Calling Non-Python External Codes
---------------------------------

Independent of its ability to process, analyze and visualize data, yt can also
serve as a mechanism for reading and selecting simulation data.  In this way,
it can be used to supply data to an external analysis routine written in
Fortran, C or C++.  This document describes how to supply that data, using the
example of a simple code that calculates the best axes that describe a
distribution of particles as a starting point.  (The underlying method is left
as an exercise for the reader; we're only currently interested in the function
specification and structs.)

If you have written a piece of code that performs some analysis function, and
you would like to include it in the base distribution of yt, we would be happy
to do so; drop us a line or see :ref:`contributing-code` for more information.

To accomplish the process of linking Python with our external code, we will be
using a language called `Cython <https://cython.org/>`_, which is
essentially a superset of Python that compiles down to C.  It is aware of NumPy
arrays, and it is able to massage data between the interpreted language Python
and C, Fortran or C++.  It will be much easier to utilize routines and analysis
code that have been separated into subroutines that accept data structures, so
we will assume that our halo axis calculator accepts a set of structs.

Our Example Code
++++++++++++++++

Here is the ``axes.h`` file in our imaginary code, which we will then wrap:

.. code-block:: c

   typedef struct structParticleCollection {
        long npart;
        double *xpos;
        double *ypos;
        double *zpos;
   } ParticleCollection;

   void calculate_axes(ParticleCollection *part,
            double *ax1, double *ax2, double *ax3);

There are several components to this analysis routine which we will have to
wrap.

#. We have to wrap the creation of an instance of ``ParticleCollection``.
#. We have to transform a set of NumPy arrays into pointers to doubles.
#. We have to create a set of doubles into which ``calculate_axes`` will be
   placing the values of the axes it calculates.
#. We have to turn the return values back into Python objects.

Each of these steps can be handled in turn, and we'll be doing it using Cython
as our interface code.

Setting Up and Building Our Wrapper
+++++++++++++++++++++++++++++++++++

To get started, we'll need to create two files:

.. code-block:: bash

   axes_calculator.pyx
   axes_calculator_setup.py

These can go anywhere, but it might be useful to put them in their own
directory.  The contents of ``axes_calculator.pyx`` will be left for the next
section, but we will need to put some boilerplate code into
``axes_calculator_setup.pyx``.  As a quick sidenote, you should call these
whatever is most appropriate for the external code you are wrapping;
``axes_calculator`` is probably not the best bet.

Here's a rough outline of what should go in ``axes_calculator_setup.py``:

.. code-block:: python

   NAME = "axes_calculator"
   EXT_SOURCES = []
   EXT_LIBRARIES = ["axes_utils", "m"]
   EXT_LIBRARY_DIRS = ["/home/rincewind/axes_calculator/"]
   EXT_INCLUDE_DIRS = []
   DEFINES = []

   from distutils.core import setup
   from distutils.extension import Extension

   from Cython.Distutils import build_ext

   ext_modules = [
       Extension(
           NAME,
           [NAME + ".pyx"] + EXT_SOURCES,
           libraries=EXT_LIBRARIES,
           library_dirs=EXT_LIBRARY_DIRS,
           include_dirs=EXT_INCLUDE_DIRS,
           define_macros=DEFINES,
       )
   ]

   setup(name=NAME, cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)

The only variables you should have to change in this are the first six, and
possibly only the first one.  We'll go through these variables one at a time.

``NAME``
   This is the name of our source file, minus the ``.pyx``.  We're also
   mandating that it be the name of the module we import.  You're free to
   modify this.
``EXT_SOURCES``
   Any additional sources can be listed here.  For instance, if you are only
   linking against a single ``.c`` file, you could list it here -- if our axes
   calculator were fully contained within a file called ``calculate_my_axes.c``
   we could link against it using this variable, and then we would not have to
   specify any libraries.  This is usually the simplest way to do things, and in
   fact, yt makes use of this itself for things like HEALPix and interpolation
   functions.
``EXT_LIBRARIES``
   Any libraries that will need to be linked against (like ``m``!) should be
   listed here.  Note that these are the name of the library minus the leading
   ``lib`` and without the trailing ``.so``.  So ``libm.so`` would become ``m``
   and ``libluggage.so`` would become ``luggage``.
``EXT_LIBRARY_DIRS``
   If the libraries listed in ``EXT_LIBRARIES`` reside in some other directory
   or directories, those directories should be listed here.  For instance,
   ``["/usr/local/lib", "/home/rincewind/luggage/"]`` .
``EXT_INCLUDE_DIRS``
   If any header files have been included that live in external directories,
   those directories should be included here.
``DEFINES``
   Any define macros that should be passed to the C compiler should be listed
   here; if they just need to be defined, then they should be specified to be
   defined as "None."  For instance, if you wanted to pass ``-DTWOFLOWER``, you
   would set this to equal: ``[("TWOFLOWER", None)]``.

To build our extension, we would run:

.. code-block:: bash

   $ python axes_calculator_setup.py build_ext -i

Note that since we don't yet have an ``axes_calculator.pyx``, this will fail.
But once we have it, it ought to run.

Writing and Calling our Wrapper
+++++++++++++++++++++++++++++++

Now we begin the tricky part, of writing our wrapper code.  We've already
figured out how to build it, which is halfway to being able to test that it
works, and we now need to start writing Cython code.

For a more detailed introduction to Cython, see the Cython documentation at
http://docs.cython.org/en/latest/ .  We'll cover a few of the basics for wrapping code
however.

To start out with, we need to open up and edit our file,
``axes_calculator.pyx``.  Open this in your favorite version of vi (mine is
vim) and we will get started by declaring the struct we need to pass in.  But
first, we need to include some header information:

.. code-block:: cython

   import numpy as np
   cimport numpy as np
   cimport cython
   from stdlib cimport malloc, free

These lines simply import and "Cython import" some common routines.  For more
information about what is already available, see the Cython documentation.  For
now, we need to start translating our data.

To do so, we tell Cython both where the struct should come from, and then we
describe the struct itself.  One fun thing to note is that if you don't need to
set or access all the values in a struct, and it just needs to be passed around
opaquely, you don't have to include them in the definition.  For an example of
this, see the ``png_writer.pyx`` file in the yt repository.  Here's the syntax
for pulling in (from a file called ``axes_calculator.h``) a struct like the one
described above:

.. code-block:: cython

   cdef extern from "axes_calculator.h":
       ctypedef struct ParticleCollection:
           long npart
           double *xpos
           double *ypos
           double *zpos

So far, pretty easy!  We've basically just translated the declaration from the
``.h`` file.  Now that we have done so, any other Cython code can create and
manipulate these ``ParticleCollection`` structs -- which we'll do shortly.
Next up, we need to declare the function we're going to call, which looks
nearly exactly like the one in the ``.h`` file.  (One common problem is that
Cython doesn't know what ``const`` means, so just remove it wherever you see
it.)  Declare it like so:

.. code-block:: cython

       void calculate_axes(ParticleCollection *part,
                double *ax1, double *ax2, double *ax3)

Note that this is indented one level, to indicate that it, too, comes from
``axes_calculator.h``.  The next step is to create a function that accepts
arrays and converts them to the format the struct likes.  We declare our
function just like we would a normal Python function, using ``def``.  You can
also use ``cdef`` if you only want to call a function from within Cython.  We
want to call it from Python, too, so we just use ``def``.  Note that we don't
here specify types for the various arguments.  In a moment we'll refine this to
have better argument types.

.. code-block:: cython

   def examine_axes(xpos, ypos, zpos):
       cdef double ax1[3], ax2[3], ax3[3]
       cdef ParticleCollection particles
       cdef int i

       particles.npart = len(xpos)
       particles.xpos = <double *> malloc(particles.npart * sizeof(double))
       particles.ypos = <double *> malloc(particles.npart * sizeof(double))
       particles.zpos = <double *> malloc(particles.npart * sizeof(double))

       for i in range(particles.npart):
           particles.xpos[i] = xpos[i]
           particles.ypos[i] = ypos[i]
           particles.zpos[i] = zpos[i]

       calculate_axes(&particles, ax1, ax2, ax3)

       free(particles.xpos)
       free(particles.ypos)
       free(particles.zpos)

       return ( (ax1[0], ax1[1], ax1[2]),
                (ax2[0], ax2[1], ax2[2]),
                (ax3[0], ax3[1], ax3[2]) )

This does the rest.  Note that we've weaved in C-type declarations (ax1, ax2,
ax3) and Python access to the variables fed in.  This function will probably be
quite slow -- because it doesn't know anything about the variables xpos, ypos,
zpos, it won't be able to speed up access to them.  Now we will see what we can
do by declaring them to be of array-type before we start handling them at all.
We can do that by annotating in the function argument list.  But first, let's
test that it works.  From the directory in which you placed these files, run:

.. code-block:: bash

   $ python2.6 setup.py build_ext -i

Now, create a sample file that feeds in the particles:

.. code-block:: python

    import axes_calculator

    axes_calculator.examine_axes(xpos, ypos, zpos)

Most of the time in that function is spent in converting the data.  So now we
can go back and we'll try again, rewriting our converter function to believe
that its being fed arrays from NumPy:

.. code-block:: cython

   def examine_axes(np.ndarray[np.float64_t, ndim=1] xpos,
                    np.ndarray[np.float64_t, ndim=1] ypos,
                    np.ndarray[np.float64_t, ndim=1] zpos):
       cdef double ax1[3], ax2[3], ax3[3]
       cdef ParticleCollection particles
       cdef int i

       particles.npart = len(xpos)
       particles.xpos = <double *> malloc(particles.npart * sizeof(double))
       particles.ypos = <double *> malloc(particles.npart * sizeof(double))
       particles.zpos = <double *> malloc(particles.npart * sizeof(double))

       for i in range(particles.npart):
           particles.xpos[i] = xpos[i]
           particles.ypos[i] = ypos[i]
           particles.zpos[i] = zpos[i]

       calculate_axes(&particles, ax1, ax2, ax3)

       free(particles.xpos)
       free(particles.ypos)
       free(particles.zpos)

       return ( (ax1[0], ax1[1], ax1[2]),
                (ax2[0], ax2[1], ax2[2]),
                (ax3[0], ax3[1], ax3[2]) )

This should be substantially faster, assuming you feed it arrays.

Now, there's one last thing we can try.  If we know our function won't modify
our arrays, and they are C-Contiguous, we can simply grab pointers to the data:

.. code-block:: cython

   def examine_axes(np.ndarray[np.float64_t, ndim=1] xpos,
                    np.ndarray[np.float64_t, ndim=1] ypos,
                    np.ndarray[np.float64_t, ndim=1] zpos):
       cdef double ax1[3], ax2[3], ax3[3]
       cdef ParticleCollection particles
       cdef int i

       particles.npart = len(xpos)
       particles.xpos = <double *> xpos.data
       particles.ypos = <double *> ypos.data
       particles.zpos = <double *> zpos.data

       for i in range(particles.npart):
           particles.xpos[i] = xpos[i]
           particles.ypos[i] = ypos[i]
           particles.zpos[i] = zpos[i]

       calculate_axes(&particles, ax1, ax2, ax3)

       return ( (ax1[0], ax1[1], ax1[2]),
                (ax2[0], ax2[1], ax2[2]),
                (ax3[0], ax3[1], ax3[2]) )

But note!  This will break or do weird things if you feed it arrays that are
non-contiguous.

At this point, you should have a mostly working piece of wrapper code.  And it
was pretty easy!  Let us know if you run into any problems, or if you are
interested in distributing your code with yt.

A complete set of files is available with this documentation.  These are
slightly different, so that the whole thing will simply compile, but they
provide a useful example.

 * `axes.c <../_static/axes.c>`_
 * `axes.h <../_static/axes.h>`_
 * `axes_calculator.pyx <../_static/axes_calculator.pyx>`_
 * `axes_calculator_setup.py <../_static/axes_calculator_setup.txt>`_

Exporting Data from yt
----------------------

yt is installed alongside h5py.  If you need to export your data from yt, to
share it with people or to use it inside another code, h5py is a good way to do
so.  You can write out complete datasets with just a few commands.  You have to
import, and then save things out into a file.

.. code-block:: python

   import h5py

   f = h5py.File("some_file.h5", mode="w")
   f.create_dataset("/data", data=some_data)

This will create ``some_file.h5`` if necessary and add a new dataset
(``/data``) to it.  Writing out in ASCII should be relatively straightforward.
For instance:

.. code-block:: python

   f = open("my_file.txt", "w")
   for halo in halos:
       x, y, z = halo.center_of_mass()
       f.write("%0.2f %0.2f %0.2f\n", x, y, z)
   f.close()

This example could be extended to work with any data object's fields, as well.

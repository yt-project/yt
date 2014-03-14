.. _parallel-computation:

Parallel Computation With YT
============================

``yt`` has been instrumented with the ability to compute many -- most, even --
quantities in parallel.  This utilizes the package 
`mpi4py <http://code.google.com/p/mpi4py>`_ to parallelize using the Message
Passing Interface, typically installed on clusters.  

.. _capabilities:

Capabilities
------------

Currently, ``yt`` is able to perform the following actions in parallel:

 * Projections (:ref:`projection-plots`)
 * Slices (:ref:`slice-plots`)
 * Cutting planes (oblique slices) (:ref:`off-axis-slices`)
 * Derived Quantities (total mass, angular momentum, etc) (:ref:`creating_derived_quantities`,
   :ref:`derived-quantities`)
 * 1-, 2-, and 3-D profiles (:ref:`generating-profiles-and-histograms`)
 * Halo finding (:ref:`halo_finding`)
 * Merger tree (:ref:`merger_tree`)
 * Two point functions (:ref:`two_point_functions`)
 * Volume rendering (:ref:`volume_rendering`)
 * Radial column density (:ref:`radial-column-density`)
 * Isocontours & flux calculations (:ref:`extracting-isocontour-information`)

This list covers just about every action ``yt`` can take!  Additionally, almost all
scripts will benefit from parallelization without any modification.  The goal
of Parallel-``yt`` has been to retain API compatibility and abstract all
parallelism.

Setting Up Parallel YT
----------------------

To run scripts in parallel, you must first install `mpi4py
<http://code.google.com/p/mpi4py>`_ as well as an MPI library, if one is not
already available on your system.  Instructions for doing so are provided on the
mpi4py website, but you may have luck by just running:

.. code-block:: bash

    $ pip install mpi4py

Once that has been installed, you're all done!  You just need to launch your 
scripts with ``mpirun`` (or equivalent) and signal to ``yt`` that you want to run 
them in parallel.  In general, that's all it takes to get a speed benefit on a 
multi-core machine.  Here is an example on an 8-core desktop:

.. code-block:: bash

    $ mpirun -np 8 python script.py --parallel

Throughout its normal operation, ``yt`` keeps you aware of what is happening with
regular messages to the stderr usually prefaced with: 

.. code-block:: bash

    yt : [INFO   ] YYY-MM-DD HH:MM:SS

However, when operating in parallel mode, ``yt`` outputs information from each
of your processors to this log mode, as in:

.. code-block:: bash

    P000 yt : [INFO   ] YYY-MM-DD HH:MM:SS
    P001 yt : [INFO   ] YYY-MM-DD HH:MM:SS

in the case of two cores being used.

It's important to note that all of the processes listed in `capabilities` work
-- and no additional work is necessary to parallelize those processes.
Furthermore, the ``yt`` command itself recognizes the ``--parallel`` option, so
those commands will work in parallel as well.

Running a ``yt`` script in parallel
-----------------------------------

Many basic ``yt`` operations will run in parallel if yt's parallelism is enabled at
startup.  For example, the following script finds the maximum density location
in the simulation and then makes a plot of the projected density:

.. code-block:: python

   from yt.pmods import *
   pf = load("RD0035/RedshiftOutput0035")
   v, c = pf.h.find_max("density")
   print v, c
   p = ProjectionPlot(pf, "x", "density")
   p.save()

If this script is run in parallel, two of the most expensive operations -
finding of the maximum density and the projection will be calulcated in
parallel.  If we save the script as ``my_script.py``, we would run it on 16 MPI
processes using the following Bash command:

.. code-block:: bash

   $ mpirun -np 16 python2.7 my_script.py --parallel

.. note::

   If you run into problems, the you can use :ref:`remote-debugging` to examine
   what went wrong.

Creating Parallel and Serial Sections in a script
+++++++++++++++++++++++++++++++++++++++++++++++++

Many ``yt`` operations will automatically run in parallel (see the next section for
a full enumeration), however some operations, particularly ones that print
output or save data to the filesystem, will be run by all processors in a
parallel script.  For example, in the script above the lines ``print v,c`` and
``p.save()`` will be run on all 16 processors.  This means that your terminal
output will contain 16 repetitions of the output of the print statement and the
plot will be saved to disk 16 times (overwritten each time).

``yt`` provides two convenience functions that make it easier to run most of a
script in parallel but run some subset of the script on only one processor.  The
first, :func:`~yt.funcs.is_root`, returns ``True`` if run on the 'root'
processor (the processor with MPI rank 0) and ``False`` otherwise.  One could
rewrite the above script to take advantage of :func:`~yt.funcs.is_root` like
so:

.. code-block:: python

   from yt.pmods import *
   pf = load("RD0035/RedshiftOutput0035")
   v, c = pf.h.find_max("density")
   p = ProjectionPlot(pf, "x", "density")
   if is_root():
       print v, c
       p.save()

The second function, :func:`~yt.funcs.only_on_root` accepts the name of a
function as well as a set of parameters and keyword arguments to pass to the
function.  This is useful when the serial component of your parallel script
would clutter the script or if you like writing your scripts as a series of
isolated function calls.  I can rewrite the example from the beginning of this
section once more using :func:`~yt.funcs.only_on_root` to give you the flavor of
how to use it:

.. code-block:: python

   from yt.pmods import *

   def print_and_save_plot(v, c, plot, print=True):
       if print:
          print v, c
       plot.save()

   pf = load("RD0035/RedshiftOutput0035")
   v, c = pf.h.find_max("density")
   p = ProjectionPlot(pf, "x", "density")
   only_on_root(print_and_save_plot, v, c, plot, print=True)

Types of Parallelism
--------------------

In order to divide up the work, ``yt`` will attempt to send different tasks to
different processors.  However, to minimize inter-process communication, YT
will decompose the information in different ways based on the task.

Spatial Decomposition
+++++++++++++++++++++

During this process, the index will be decomposed along either all three
axes or along an image plane, if the process is that of projection.  This type
of parallelism is overall less efficient than grid-based parallelism, but it
has been shown to obtain good results overall.

The following operations use spatial decomposition:

  * Halo finding
  * Merger tree
  * Two point functions
  * Volume rendering
  * Radial column density

Grid Decomposition
++++++++++++++++++

The alternative to spatial decomposition is a simple round-robin of the grids.
This process allows ``yt`` to pool data access to a given Enzo data file, which
ultimately results in faster read times and better parallelism.

The following operations use grid decomposition:

  * Projections
  * Slices
  * Cutting planes
  * Derived Quantities
  * 1-, 2-, and 3-D profiles
  * Isocontours & flux calculations

Object-Based
++++++++++++

In a fashion similar to grid decomposition, computation can be parallelized
over objects. This is especially useful for
`embarrassingly parallel <http://en.wikipedia.org/wiki/Embarrassingly_parallel>`_
tasks where the items to be worked on can be split into separate chunks and
saved to a list. The list is then split up and each MPI task performs parts of
it independently.

.. _parallelizing-your-analysis:

Parallelizing Your Analysis
---------------------------

It is easy within ``yt`` to parallelize a list of tasks, as long as those tasks
are independent of one another.
Using object-based parallelism, the function :func:`parallel_objects` will
automatically split up a list of tasks over the specified number of processors
(or cores).
Please see this heavily-commented example:

.. code-block:: python
   
   # As always...
   from yt.mods import *
   
   import glob
   
   # The number 4, below, is the number of processes to parallelize over, which
   # is generally equal to the number of MPI tasks the job is launched with.
   # If num_procs is set to zero or a negative number, the for loop below
   # will be run such that each iteration of the loop is done by a single MPI
   # task. Put another way, setting it to zero means that no matter how many
   # MPI tasks the job is run with, num_procs will default to the number of
   # MPI tasks automatically.
   num_procs = 4
   
   # fns is a list of all the simulation data files in the current directory.
   fns = glob.glob("./plot*")
   fns.sort()

   # This dict will store information collected in the loop, below.
   # Inside the loop each task will have a local copy of the dict, but
   # the dict will be combined once the loop finishes.
   my_storage = {}

   # In this example, because the storage option is used in the
   # parallel_objects function, the loop yields a tuple, which gets used
   # as (sto, fn) inside the loop.
   # In the loop, sto is essentially my_storage, but a local copy of it.
   # If data does not need to be combined after the loop is done, the line
   # would look like:
   #       for fn in parallel_objects(fns, num_procs):
   for sto, fn in parallel_objects(fns, num_procs, storage = my_storage):

       # Open a data file, remembering that fn is different on each task.
       pf = load(fn)
       dd = pf.h.all_data()

       # This copies fn and the min/max of density to the local copy of
       # my_storage
       sto.result_id = fn
       sto.result = dd.quantities["Extrema"]("density")

       # Makes and saves a plot of the gas density.
       p = ProjectionPlot(pf, "x", "density")
       p.save()

   # At this point, as the loop exits, the local copies of my_storage are
   # combined such that all tasks now have an identical and full version of
   # my_storage. Until this point, each task is unaware of what the other
   # tasks have produced.
   # Below, the values in my_storage are printed by only one task. The other
   # tasks do nothing.
   if is_root()
       for fn, vals in sorted(my_storage.items()):
           print fn, vals

This example above can be modified to loop over anything that can be saved to
a Python list: halos, data files, arrays, and more.

.. _parallel-time-series-analysis:

Parallel Time Series Analysis
-----------------------------

The same :func:`parallel_objects` machinery discussed above is turned on by
default when using a ``DatasetSeries`` object (see :ref:`time-series-analysis`)
to iterate over simulation outputs.  The syntax for this is very simple.  As an
example, we can use the following script to find the angular momentum vector in
a 1 pc sphere centered on the maximum density cell in a large number of
simulation outputs:

.. code-block:: python

   from yt.pmods import *
   ts = DatasetSeries.from_filenames("DD*/output_*", parallel = True)
   sphere = ts.sphere("max", (1.0, "pc"))
   L_vecs = sphere.quantities["AngularMomentumVector"]()

Note that this script can be run in serial or parallel with an arbitrary number
of processors.  When running in parallel, each output is given to a different
processor.  By default, parallel is set to ``True``, so you do not have to
explicitly set ``parallel = True`` as in the above example. 

One could get the same effect by iterating over the individual parameter files
in the DatasetSeries object:

.. code-block:: python

   from yt.pmods import *
   ts = DatasetSeries.from_filenames("DD*/output_*", parallel = True)
   my_storage = {}
   for sto,pf in ts.piter(storage=my_storage):
       sphere = pf.sphere("max", (1.0, "pc"))
       L_vec = sphere.quantities["AngularMomentumVector"]()
       sto.result_id = pf.parameter_filename
       sto.result = L_vec

   L_vecs = []
   for fn, L_vec in sorted(my_storage.items()):
       L_vecs.append(L_vec)


You can also request a fixed number of processors to calculate each
angular momentum vector.  For example, this script will calculate each angular
momentum vector using 4 workgroups, splitting up the pool available processors.
Note that parallel=1 implies that the analysis will be run using 1 workgroup, 
whereas parallel=True will run with Nprocs workgroups.

.. code-block:: python

   from yt.pmods import *
   ts = DatasetSeries.from_filenames("DD*/output_*", parallel = 4)
   sphere = ts.sphere("max", (1.0, "pc))
   L_vecs = sphere.quantities["AngularMomentumVector"]()

If you do not want to use ``parallel_objects`` parallelism when using a
TimeSeries object, set ``parallel = False``.  When running python in parallel,
this will use all of the available processors to evaluate the requested
operation on each simulation output.  Some care and possibly trial and error
might be necessary to estimate the correct settings for your Simulation
outputs.

Parallel Performance, Resources, and Tuning
-------------------------------------------

Optimizing parallel jobs in ``yt`` is difficult; there are many parameters that
affect how well and quickly the job runs.  In many cases, the only way to find
out what the minimum (or optimal) number of processors is, or amount of memory
needed, is through trial and error.  However, this section will attempt to
provide some insight into what are good starting values for a given parallel
task.

Grid Decomposition
++++++++++++++++++

In general, these types of parallel calculations scale very well with number of
processors.
They are also fairly memory-conservative.
The two limiting factors is therefore the number of grids in the dataset,
and the speed of the disk the data is stored on.
There is no point in running a parallel job of this kind with more processors
than grids, because the extra processors will do absolutely nothing, and will
in fact probably just serve to slow down the whole calculation due to the extra
overhead.
The speed of the disk is also a consideration - if it is not a high-end parallel
file system, adding more tasks will not speed up the calculation if the disk
is already swamped with activity.

The best advice for these sort of calculations is to run 
with just a few processors and go from there, seeing if it the runtime
improves noticeably.

**Projections, Slices, and Cutting Planes**

Projections, slices and cutting planes are the most common methods of creating
two-dimensional representations of data.  All three have been parallelized in a
grid-based fashion.

 * Projections: projections are parallelized utilizing a quad-tree approach.
   Data is loaded for each processor, typically by a process that consolidates
   open/close/read operations, and each grid is then iterated over and cells
   are deposited into a data structure that stores values corresponding to
   positions in the two-dimensional plane.  This provides excellent load
   balancing, and in serial is quite fast.  However, as of ``yt`` 2.3, the
   operation by which quadtrees are joined across processors scales poorly;
   while memory consumption scales well, the time to completion does not.  As
   such, projections can often be done very fast when operating only on a single
   processor!  The quadtree algorithm can be used inline (and, indeed, it is
   for this reason that it is slow.)  It is recommended that you attempt to
   project in serial before projecting in parallel; even for the very largest
   datasets (Enzo 1024^3 root grid with 7 levels of refinement) in the absence
   of IO the quadtree algorithm takes only three minutes or so on a decent
   processor.
 * Slices: to generate a slice, grids that intersect a given slice are iterated
   over and their finest-resolution cells are deposited.  The grids are
   decomposed via standard load balancing.  While this operation is parallel,
   **it is almost never necessary to slice a dataset in parallel**, as all data is
   loaded on demand anyway.  The slice operation has been parallelized so as to
   enable slicing when running *in situ*.
 * Cutting planes: cutting planes are parallelized exactly as slices are.
   However, in contrast to slices, because the data-selection operation can be
   much more time consuming, cutting planes often benefit from parallelism.

Object-Based
++++++++++++

Like grid decomposition, it does not help to run with more processors than the
number of objects to be iterated over.
There is also the matter of the kind of work being done on each object, and
whether it is disk-intensive, cpu-intensive, or memory-intensive.
It is up to the user to figure out what limits the performance of their script,
and use the correct amount of resources, accordingly.

Disk-intensive jobs are limited by the speed of the file system, as above,
and extra processors beyond its capability are likely counter-productive.
It may require some testing or research (e.g. supercomputer documentation)
to find out what the file system is capable of.

If it is cpu-intensive, it's best to use as many processors as possible
and practical.

For a memory-intensive job, each processor needs to be able to allocate enough
memory, which may mean using fewer than the maximum number of tasks per compute
node, and increasing the number of nodes.
The memory used per processor should be calculated, compared to the memory
on each compute node, which dictates how many tasks per node.
After that, the number of processors used overall is dictated by the 
disk system or CPU-intensity of the job.


Domain Decomposition
++++++++++++++++++++

The various types of analysis that utilize domain decomposition use them in
different enough ways that they are be discussed separately.

**Halo-Finding**

Halo finding, along with the merger tree that uses halo finding, operates
on the particles in the volume, and is therefore mostly grid-agnostic.
Generally, the biggest concern for halo finding is the amount of memory needed.
There is subtle art in estimating the amount of memory needed for halo finding,
but a rule of thumb is that Parallel HOP (:func:`parallelHF`) is the most
memory-intensive, followed by plain HOP (:func:`HaloFinder`),
with Friends of Friends (:func:`FOFHaloFinder`) being
the most memory-conservative.
It has been found that :func:`parallelHF` needs roughly
1 MB of memory per 5,000
particles, although recent work has improved this and the memory requirement
is now smaller than this. But this is a good starting point for beginning to
calculate the memory required for halo-finding.

**Two point functions**

Please see :ref:`tpf_strategies` for more details.

**Volume Rendering**

The simplest way to think about volume rendering, and the radial column density
module that uses it, is that it load-balances over the grids in the dataset.
Each processor is given roughly the same sized volume to operate on.
In practice, there are just a few things to keep in mind when doing volume
rendering.
First, it only uses a power of two number of processors.
If the job is run with 100 processors, only 64 of them will actually do anything.
Second, the absolute maximum number of processors is the number of grids.
But in order to keep work distributed evenly, typically the number of processors
should be no greater than one-eighth or one-quarter the number of processors
that were used to produce the dataset.

Additional Tips
---------------

  * Don't be afraid to change how a parallel job is run. Change the
    number of processors, or memory allocated, and see if things work better
    or worse. After all, it's just a computer, it doesn't pass moral judgment!

  * Similarly, human time is more valuable than computer time. Try increasing
    the number of processors, and see if the runtime drops significantly.
    There will be a sweet spot between speed of run and the waiting time in
    the job scheduler queue; it may be worth trying to find it.

  * If you are using object-based parallelism but doing CPU-intensive computations
    on each object, you may find that setting ``num_procs`` equal to the 
    number of processors per compute node can lead to significant speedups.
    By default, most mpi implementations will assign tasks to processors on a
    'by-slot' basis, so this setting will tell ``yt`` to do computations on a single
    object using only the processors on a single compute node.  A nice application
    for this type of parallelism is calculating a list of derived quantities for 
    a large number of simulation outputs.

  * It is impossible to tune a parallel operation without understanding what's
    going on. Read the documentation, look at the underlying code, or talk to
    other ``yt`` users. Get informed!
    
  * Sometimes it is difficult to know if a job is cpu, memory, or disk
    intensive, especially if the parallel job utilizes several of the kinds of
    parallelism discussed above. In this case, it may be worthwhile to put
    some simple timers in your script (as below) around different parts.
    
    .. code-block:: python
    
       from yt.mods import *
       import time
       
       pf = load("DD0152")
       t0 = time.time()
       bigstuff, hugestuff = StuffFinder(pf)
       BigHugeStuffParallelFunction(pf, bigstuff, hugestuff)
       t1 = time.time()
       for i in range(1000000):
           tinystuff, ministuff = GetTinyMiniStuffOffDisk("in%06d.txt" % i)
           array = TinyTeensyParallelFunction(pf, tinystuff, ministuff)
           SaveTinyMiniStuffToDisk("out%06d.txt" % i, array)
       t2 = time.time()
       
       if is_root()
           print "BigStuff took %.5e sec, TinyStuff took %.5e sec" % (t1 - t0, t2 - t1)
  
  * Remember that if the script handles disk IO explicitly, and does not use
    a built-in ``yt`` function to write data to disk,
    care must be taken to
    avoid `race-conditions <http://en.wikipedia.org/wiki/Race_conditions>`_.
    Be explicit about which MPI task writes to disk using a construction
    something like this:
    
    .. code-block:: python
       
       if is_root()
           file = open("out.txt", "w")
           file.write(stuff)
           file.close()

  * Many supercomputers allow users to ssh into the nodes that their job is
    running on.
    Many job schedulers send the names of the nodes that are
    used in the notification emails, or a command like ``qstat -f NNNN``, where
    ``NNNN`` is the job ID, will also show this information.
    By ssh-ing into nodes, the memory usage of each task can be viewed in
    real-time as the job runs (using ``top``, for example),
    and can give valuable feedback about the
    resources the task requires.
    
An Advanced Worked Example
--------------------------

Below is a script used to calculate the redshift of first 99.9% ionization in a
simulation.  This script was designed to analyze a set of 100 outputs on
Gordon, running on 128 processors.  This script goes through three phases:

 #. Define a new derived field, which calculates the fraction of ionized
    hydrogen as a function only of the total hydrogen density.
 #. Load a time series up, specifying ``parallel = 8``.  This means that it
    will decompose into 8 jobs.  So if we ran on 128 processors, we would have
    16 processors assigned to each output in the time series.
 #. Creating a big cube that will hold our results for this set of processors.
    Note that this will be only for each output considered by this processor,
    and this cube will not necessarily be filled in in every cell.
 #. For each output, distribute the grids to each of the sixteen processors
    working on that output.  Each of these takes the max of the ionized
    redshift in their zone versus the accumulation cube.
 #. Iterate over slabs and find the maximum redshift in each slab of our
    accumulation cube.

At the end, the root processor (of the global calculation) writes out an
ionization cube that contains the redshift of first reionization for each zone
across all outputs.

.. literalinclude:: ionization_cube.py

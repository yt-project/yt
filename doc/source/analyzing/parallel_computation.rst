.. _parallel-computation:

Parallel Computation With yt
============================

yt has been instrumented with the ability to compute many -- most, even --
quantities in parallel.  This utilizes the package
`mpi4py <https://bitbucket.org/mpi4py/mpi4py>`_ to parallelize using the Message
Passing Interface, typically installed on clusters.

.. _capabilities:

Capabilities
------------

Currently, yt is able to perform the following actions in parallel:

* Projections (:ref:`projection-plots`)
* Slices (:ref:`slice-plots`)
* Cutting planes (oblique slices) (:ref:`off-axis-slices`)
* Covering grids (:ref:`examining-grid-data-in-a-fixed-resolution-array`)
* Derived Quantities (total mass, angular momentum, etc) (:ref:`creating_derived_quantities`,
  :ref:`derived-quantities`)
* 1-, 2-, and 3-D profiles (:ref:`generating-profiles-and-histograms`)
* Halo analysis (:ref:`halo-analysis`)
* Volume rendering (:ref:`volume_rendering`)
* Isocontours & flux calculations (:ref:`extracting-isocontour-information`)

This list covers just about every action yt can take!  Additionally, almost all
scripts will benefit from parallelization with minimal modification.  The goal
of Parallel-yt has been to retain API compatibility and abstract all
parallelism.

Setting Up Parallel yt
--------------------------

To run scripts in parallel, you must first install `mpi4py
<https://bitbucket.org/mpi4py/mpi4py>`_ as well as an MPI library, if one is not
already available on your system.  Instructions for doing so are provided on the
mpi4py website, but you may have luck by just running:

.. code-block:: bash

    $ pip install mpi4py

If you have an Anaconda installation of yt and there is no MPI library on the
system you are using try:

.. code-block:: bash

    $ conda install mpi4py

This will install `MPICH2 <https://www.mpich.org/>`_ and will interfere with
other MPI libraries that are already installed. Therefore, it is preferable to
use the ``pip`` installation method.

Once mpi4py has been installed, you're all done!  You just need to launch your
scripts with ``mpirun`` (or equivalent) and signal to yt that you want to
run them in parallel by invoking the ``yt.enable_parallelism()`` function in
your script.  In general, that's all it takes to get a speed benefit on a
multi-core machine.  Here is an example on an 8-core desktop:

.. code-block:: bash

    $ mpirun -np 8 python script.py

Throughout its normal operation, yt keeps you aware of what is happening with
regular messages to the stderr usually prefaced with:

.. code-block:: bash

    yt : [INFO   ] YYY-MM-DD HH:MM:SS

However, when operating in parallel mode, yt outputs information from each
of your processors to this log mode, as in:

.. code-block:: bash

    P000 yt : [INFO   ] YYY-MM-DD HH:MM:SS
    P001 yt : [INFO   ] YYY-MM-DD HH:MM:SS

in the case of two cores being used.

It's important to note that all of the processes listed in :ref:`capabilities`
work in parallel -- and no additional work is necessary to parallelize those
processes.

Running a yt Script in Parallel
-------------------------------

Many basic yt operations will run in parallel if yt's parallelism is enabled at
startup.  For example, the following script finds the maximum density location
in the simulation and then makes a plot of the projected density:

.. code-block:: python

   import yt

   yt.enable_parallelism()

   ds = yt.load("RD0035/RedshiftOutput0035")
   v, c = ds.find_max("density")
   print(v, c)
   p = yt.ProjectionPlot(ds, "x", "density")
   p.save()

If this script is run in parallel, two of the most expensive operations -
finding of the maximum density and the projection will be calculated in
parallel.  If we save the script as ``my_script.py``, we would run it on 16 MPI
processes using the following Bash command:

.. code-block:: bash

   $ mpirun -np 16 python my_script.py

.. note::

   If you run into problems, the you can use :ref:`remote-debugging` to examine
   what went wrong.

How do I run my yt job on a subset of available processes
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can set the ``communicator`` keyword in the
:func:`~yt.utilities.parallel_tools.parallel_analysis_interface.enable_parallelism`
call to a specific MPI communicator to specify a subset of available MPI
processes.  If none is specified, it defaults to ``COMM_WORLD``.

Creating Parallel and Serial Sections in a Script
+++++++++++++++++++++++++++++++++++++++++++++++++

Many yt operations will automatically run in parallel (see the next section for
a full enumeration), however some operations, particularly ones that print
output or save data to the filesystem, will be run by all processors in a
parallel script.  For example, in the script above the lines ``print(v, c)`` and
``p.save()`` will be run on all 16 processors.  This means that your terminal
output will contain 16 repetitions of the output of the print statement and the
plot will be saved to disk 16 times (overwritten each time).

yt provides two convenience functions that make it easier to run most of a
script in parallel but run some subset of the script on only one processor.  The
first, :func:`~yt.funcs.is_root`, returns ``True`` if run on the 'root'
processor (the processor with MPI rank 0) and ``False`` otherwise.  One could
rewrite the above script to take advantage of :func:`~yt.funcs.is_root` like
so:

.. code-block:: python

   import yt

   yt.enable_parallelism()

   ds = yt.load("RD0035/RedshiftOutput0035")
   v, c = ds.find_max("density")
   p = yt.ProjectionPlot(ds, "x", "density")
   if yt.is_root():
       print(v, c)
       p.save()

The second function, :func:`~yt.funcs.only_on_root` accepts the name of a
function as well as a set of parameters and keyword arguments to pass to the
function.  This is useful when the serial component of your parallel script
would clutter the script or if you like writing your scripts as a series of
isolated function calls.  I can rewrite the example from the beginning of this
section once more using :func:`~yt.funcs.only_on_root` to give you the flavor of
how to use it:

.. code-block:: python

   import yt

   yt.enable_parallelism()


   def print_and_save_plot(v, c, plot, verbose=True):
       if verbose:
           print(v, c)
       plot.save()


   ds = yt.load("RD0035/RedshiftOutput0035")
   v, c = ds.find_max("density")
   p = yt.ProjectionPlot(ds, "x", "density")
   yt.only_on_root(print_and_save_plot, v, c, plot, verbose=True)

Types of Parallelism
--------------------

In order to divide up the work, yt will attempt to send different tasks to
different processors.  However, to minimize inter-process communication, yt
will decompose the information in different ways based on the task.

Spatial Decomposition
+++++++++++++++++++++

During this process, the index will be decomposed along either all three
axes or along an image plane, if the process is that of projection.  This type
of parallelism is overall less efficient than grid-based parallelism, but it
has been shown to obtain good results overall.

The following operations use spatial decomposition:

* :ref:`halo-analysis`
* :ref:`volume_rendering`

Grid Decomposition
++++++++++++++++++

The alternative to spatial decomposition is a simple round-robin of data chunks,
which could be grids, octs, or whatever chunking mechanism is used by the code
frontend begin used.  This process allows yt to pool data access to a given
data file, which ultimately results in faster read times and better parallelism.

The following operations use chunk decomposition:

* Projections (see :ref:`available-objects`)
* Slices (see :ref:`available-objects`)
* Cutting planes (see :ref:`available-objects`)
* Covering grids (see :ref:`construction-objects`)
* Derived Quantities (see :ref:`derived-quantities`)
* 1-, 2-, and 3-D profiles (see :ref:`generating-profiles-and-histograms`)
* Isocontours & flux calculations (see :ref:`surfaces`)

Parallelization over Multiple Objects and Datasets
++++++++++++++++++++++++++++++++++++++++++++++++++

If you have a set of computational steps that need to apply identically and
independently to several different objects or datasets, a so-called
`embarrassingly parallel <https://en.wikipedia.org/wiki/Embarrassingly_parallel>`_
task, yt can do that easily.  See the sections below on
:ref:`parallelizing-your-analysis` and :ref:`parallel-time-series-analysis`.

Use of ``piter()``
^^^^^^^^^^^^^^^^^^

If you use parallelism over objects or datasets, you will encounter
the :func:`~yt.data_objects.time_series.DatasetSeries.piter` function.
:func:`~yt.data_objects.time_series.DatasetSeries.piter` is a parallel iterator,
which effectively doles out each item of a DatasetSeries object to a different
processor.  In serial processing, you might iterate over a DatasetSeries by:

.. code-block:: python

    for dataset in dataset_series:
        ...  # process

But in parallel, you can use ``piter()`` to force each dataset to go to
a different processor:

.. code-block:: python

    yt.enable_parallelism()
    for dataset in dataset_series.piter():
        ...  # process

In order to store information from the parallel processing step to
a data structure that exists on all of the processors operating in parallel
we offer the ``storage`` keyword in the
:func:`~yt.data_objects.time_series.DatasetSeries.piter` function.
You may define an empty dictionary and include it as the keyword argument
``storage`` to :func:`~yt.data_objects.time_series.DatasetSeries.piter`.
Then, during the processing step, you can access
this dictionary as the ``sto`` object.  After the
loop is finished, the dictionary is re-aggregated from all of the processors,
and you can access the contents:

.. code-block:: python

    yt.enable_parallelism()
    my_dictionary = {}
    for sto, dataset in dataset_series.piter(storage=my_dictionary):
        ...  # process
        sto.result = ...  # some information processed for this dataset
        sto.result_id = ...  # some identifier for this dataset

    print(my_dictionary)

By default, the dataset series will be divided as equally as possible
among the cores.  Often some datasets will require more work than
others.  We offer the ``dynamic`` keyword in the
:func:`~yt.data_objects.time_series.DatasetSeries.piter` function to
enable dynamic load balancing with a task queue.  Dynamic load
balancing works best with more cores and a variable workload.  Here
one process will act as a server to assign the next available dataset
to any free client.  For example, a 16 core job will have 15 cores
analyzing the data with 1 core acting as the task manager.

.. _parallelizing-your-analysis:

Parallelizing over Multiple Objects
-----------------------------------

It is easy within yt to parallelize a list of tasks, as long as those tasks
are independent of one another. Using object-based parallelism, the function
:func:`~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_objects`
will automatically split up a list of tasks over the specified number of
processors (or cores).  Please see this heavily-commented example:

.. code-block:: python

   # As always...
   import yt

   yt.enable_parallelism()

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
   for sto, fn in yt.parallel_objects(fns, num_procs, storage=my_storage):

       # Open a data file, remembering that fn is different on each task.
       ds = yt.load(fn)
       dd = ds.all_data()

       # This copies fn and the min/max of density to the local copy of
       # my_storage
       sto.result_id = fn
       sto.result = dd.quantities.extrema("density")

       # Makes and saves a plot of the gas density.
       p = yt.ProjectionPlot(ds, "x", "density")
       p.save()

   # At this point, as the loop exits, the local copies of my_storage are
   # combined such that all tasks now have an identical and full version of
   # my_storage. Until this point, each task is unaware of what the other
   # tasks have produced.
   # Below, the values in my_storage are printed by only one task. The other
   # tasks do nothing.
   if yt.is_root():
       for fn, vals in sorted(my_storage.items()):
           print(fn, vals)

This example above can be modified to loop over anything that can be saved to
a Python list: halos, data files, arrays, and more.

.. _parallel-time-series-analysis:

Parallelization over Multiple Datasets (including Time Series)
--------------------------------------------------------------

The same ``parallel_objects`` machinery discussed above is turned on by
default when using a :class:`~yt.data_objects.time_series.DatasetSeries` object
(see :ref:`time-series-analysis`) to iterate over simulation outputs.  The
syntax for this is very simple.  As an example, we can use the following script
to find the angular momentum vector in a 1 pc sphere centered on the maximum
density cell in a large number of simulation outputs:

.. code-block:: python

   import yt

   yt.enable_parallelism()

   # Load all of the DD*/output_* files into a DatasetSeries object
   # in this case it is a Time Series
   ts = yt.load("DD*/output_*")

   # Define an empty storage dictionary for collecting information
   # in parallel through processing
   storage = {}

   # Use piter() to iterate over the time series, one proc per dataset
   # and store the resulting information from each dataset in
   # the storage dictionary
   for sto, ds in ts.piter(storage=storage):
       sphere = ds.sphere("max", (1.0, "pc"))
       sto.result = sphere.quantities.angular_momentum_vector()
       sto.result_id = str(ds)

   # Print out the angular momentum vector for all of the datasets
   for L in sorted(storage.items()):
       print(L)

Note that this script can be run in serial or parallel with an arbitrary number
of processors.  When running in parallel, each output is given to a different
processor.

You can also request a fixed number of processors to calculate each
angular momentum vector.  For example, the following script will calculate each
angular momentum vector using 4 workgroups, splitting up the pool available
processors.  Note that parallel=1 implies that the analysis will be run using
1 workgroup, whereas parallel=True will run with Nprocs workgroups.

.. code-block:: python

   import yt

   yt.enable_parallelism()

   ts = yt.DatasetSeries("DD*/output_*", parallel=4)

   for ds in ts.piter():
       sphere = ds.sphere("max", (1.0, "pc"))
       L_vecs = sphere.quantities.angular_momentum_vector()

If you do not want to use ``parallel_objects`` parallelism when using a
DatasetSeries object, set ``parallel = False``.  When running python in parallel,
this will use all of the available processors to evaluate the requested
operation on each simulation output.  Some care and possibly trial and error
might be necessary to estimate the correct settings for your simulation
outputs.

Note, when iterating over several large datasets, running out of memory may
become an issue as the internal data structures associated with each dataset
may not be properly de-allocated at the end of an iteration. If memory use
becomes a problem, it may be necessary to manually delete some of the larger
data structures.

.. code-block:: python

   import yt

   yt.enable_parallelism()

   ts = yt.DatasetSeries("DD*/output_*", parallel=4)

   for ds in ts.piter():
       # do analysis here

       ds.index.clear_all_data()

Multi-level Parallelism
-----------------------

By default, the
:func:`~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_objects`
and :func:`~yt.data_objects.time_series.DatasetSeries.piter` functions will allocate a
single processor to each iteration of the parallelized loop. However, there may be
situations in which it is advantageous to have multiple processors working together
on each loop iteration. Like with any traditional for loop, nested loops with multiple
calls to :func:`~yt.utilities.parallel_tools.parallel_analysis_interface.enable_parallelism`
can be used to parallelize the functionality within a given loop iteration.

In the example below, we will create projections along the x, y, and z axis of the
density and temperature fields. We will assume a total of 6 processors are available,
allowing us to allocate to processors to each axis and project each field with a
separate processor.

.. code-block:: python

   import yt

   yt.enable_parallelism()

   # assume 6 total cores
   # allocate 3 work groups of 2 cores each
   for ax in yt.parallel_objects("xyz", njobs=3):

       # project each field with one of the two cores in the workgroup
       for field in yt.parallel_objects(["density", "temperature"]):
           p = yt.ProjectionPlot(ds, ax, field, weight_field="density")
           p.save("figures/")

Note, in the above example, if the inner
:func:`~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_objects`
call were removed from the loop, the two-processor work group would work together to
project each of the density and temperature fields. This is because the projection
functionality itself is parallelized internally.

The :func:`~yt.data_objects.time_series.DatasetSeries.piter` function can also be used
in the above manner with nested
:func:`~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_objects`
loops to allocate multiple processors to work on each dataset. As discussed above in
:ref:`parallel-time-series-analysis`, the ``parallel`` keyword is used to control
the number of workgroups created for iterating over multiple datasets.

Parallel Performance, Resources, and Tuning
-------------------------------------------

Optimizing parallel jobs in yt is difficult; there are many parameters that
affect how well and quickly the job runs.  In many cases, the only way to find
out what the minimum (or optimal) number of processors is, or amount of memory
needed, is through trial and error.  However, this section will attempt to
provide some insight into what are good starting values for a given parallel
task.

Chunk Decomposition
+++++++++++++++++++

In general, these types of parallel calculations scale very well with number of
processors.  They are also fairly memory-conservative.  The two limiting factors
is therefore the number of chunks in the dataset, and the speed of the disk the
data is stored on.  There is no point in running a parallel job of this kind
with more processors than chunks, because the extra processors will do absolutely
nothing, and will in fact probably just serve to slow down the whole calculation
due to the extra overhead.  The speed of the disk is also a consideration - if
it is not a high-end parallel file system, adding more tasks will not speed up
the calculation if the disk is already swamped with activity.

The best advice for these sort of calculations is to run with just a few
processors and go from there, seeing if it the runtime improves noticeably.

**Projections, Slices, Cutting Planes and Covering Grids**

Projections, slices and cutting planes are the most common methods of creating
two-dimensional representations of data.  All three have been parallelized in a
chunk-based fashion.

* **Projections**: projections are parallelized utilizing a quad-tree approach.
  Data is loaded for each processor, typically by a process that consolidates
  open/close/read operations, and each grid is then iterated over and cells are
  deposited into a data structure that stores values corresponding to positions
  in the two-dimensional plane.  This provides excellent load balancing, and in
  serial is quite fast.  However, the operation by which quadtrees are joined
  across processors scales poorly; while memory consumption scales well, the
  time to completion does not.  As such, projections can often be done very
  fast when operating only on a single processor!  The quadtree algorithm can
  be used inline (and, indeed, it is for this reason that it is slow.)  It is
  recommended that you attempt to project in serial before projecting in
  parallel; even for the very largest datasets (Enzo 1024^3 root grid with 7
  levels of refinement) in the absence of IO the quadtree algorithm takes only
  three minutes or so on a decent processor.

* **Slices**: to generate a slice, chunks that intersect a given slice are iterated
  over and their finest-resolution cells are deposited.  The chunks are
  decomposed via standard load balancing.  While this operation is parallel,
  **it is almost never necessary to slice a dataset in parallel**, as all data is
  loaded on demand anyway.  The slice operation has been parallelized so as to
  enable slicing when running *in situ*.

* **Cutting planes**: cutting planes are parallelized exactly as slices are.
  However, in contrast to slices, because the data-selection operation can be
  much more time consuming, cutting planes often benefit from parallelism.

* **Covering Grids**: covering grids are parallelized exactly as slices are.

Object-Based
++++++++++++

Like chunk decomposition, it does not help to run with more processors than the
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
different enough ways that they are discussed separately.

**Halo-Finding**

Halo finding, along with the merger tree that uses halo finding, operates on the
particles in the volume, and is therefore mostly chunk-agnostic.  Generally, the
biggest concern for halo finding is the amount of memory needed.  There is
subtle art in estimating the amount of memory needed for halo finding, but a
rule of thumb is that the HOP halo finder is the most memory intensive
(:func:`HaloFinder`), and Friends of Friends (:func:`FOFHaloFinder`) being the
most memory-conservative. For more information, see :ref:`halo-analysis`.

**Volume Rendering**

The simplest way to think about volume rendering, is that it load-balances over
the i/o chunks in the dataset.  Each processor is given roughly the same sized
volume to operate on.  In practice, there are just a few things to keep in mind
when doing volume rendering.  First, it only uses a power of two number of
processors.  If the job is run with 100 processors, only 64 of them will
actually do anything.  Second, the absolute maximum number of processors is the
number of chunks.  In order to keep work distributed evenly, typically the
number of processors should be no greater than one-eighth or one-quarter the
number of processors that were used to produce the dataset.
For more information, see :ref:`volume_rendering`.

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
  'by-slot' basis, so this setting will tell yt to do computations on a single
  object using only the processors on a single compute node.  A nice application
  for this type of parallelism is calculating a list of derived quantities for
  a large number of simulation outputs.

* It is impossible to tune a parallel operation without understanding what's
  going on. Read the documentation, look at the underlying code, or talk to
  other yt users. Get informed!

* Sometimes it is difficult to know if a job is cpu, memory, or disk
  intensive, especially if the parallel job utilizes several of the kinds of
  parallelism discussed above. In this case, it may be worthwhile to put
  some simple timers in your script (as below) around different parts.

.. code-block:: python

   import time

   import yt

   yt.enable_parallelism()

   ds = yt.load("DD0152")
   t0 = time.time()
   bigstuff, hugestuff = StuffFinder(ds)
   BigHugeStuffParallelFunction(ds, bigstuff, hugestuff)
   t1 = time.time()
   for i in range(1000000):
       tinystuff, ministuff = GetTinyMiniStuffOffDisk("in%06d.txt" % i)
       array = TinyTeensyParallelFunction(ds, tinystuff, ministuff)
       SaveTinyMiniStuffToDisk("out%06d.txt" % i, array)
   t2 = time.time()

   if yt.is_root():
       print(
           "BigStuff took {:.5e} sec, TinyStuff took {:.5e} sec".format(t1 - t0, t2 - t1)
       )

* Remember that if the script handles disk IO explicitly, and does not use
  a built-in yt function to write data to disk,
  care must be taken to
  avoid `race-conditions <https://en.wikipedia.org/wiki/Race_conditions>`_.
  Be explicit about which MPI task writes to disk using a construction
  something like this:

.. code-block:: python

   if yt.is_root():
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
   and this cube will not necessarily be filled in every cell.
#. For each output, distribute the grids to each of the sixteen processors
   working on that output.  Each of these takes the max of the ionized
   redshift in their zone versus the accumulation cube.
#. Iterate over slabs and find the maximum redshift in each slab of our
   accumulation cube.

At the end, the root processor (of the global calculation) writes out an
ionization cube that contains the redshift of first reionization for each zone
across all outputs.

.. literalinclude:: ionization_cube.py

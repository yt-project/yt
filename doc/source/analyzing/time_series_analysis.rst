.. _time-series-analysis:

Time Series Analysis
====================

Often, one wants to analyze a continuous set of outputs from a simulation in a
uniform manner.  A simple example would be to calculate the peak density in a
set of outputs that were written out.  The problem with time series analysis in
yt is general an issue of verbosity and clunkiness. Typically, one sets up a 
loop:

.. code-block:: python

   for dsi in range(30):
       fn = "DD%04i/DD%04i" % (dsi, dsi)
       ds = load(fn)
       process_output(ds)

But this is not really very nice.  This ends up requiring a lot of maintenance.
The :class:`~yt.data_objects.time_series.DatasetSeries` object has been
designed to remove some of this clunkiness and present an easier, more unified
approach to analyzing sets of data.  Even better,
:class:`~yt.data_objects.time_series.DatasetSeries` works in parallel by
default (see :ref:`parallel-computation`), so you can use a ``DatasetSeries``
object to quickly and easily parallelize your analysis.  Since doing the same
analysis task on many simulation outputs is 'embarrassingly' parallel, this
naturally allows for almost arbitrary speedup - limited only by the number of
available processors and the number of simulation outputs.

The idea behind the current implementation of time series analysis is that
the underlying data and the operators that act on that data can and should be
distinct.  There are several operators provided, as well as facilities for
creating your own, and these operators can be applied either to datasets on the
whole or to subregions of individual datasets.

The simplest mechanism for creating a ``DatasetSeries`` object is to use the
class method
:meth:`~yt.data_objects.time_series.DatasetSeries.from_filenames`.  This
method accepts a list of strings that can be supplied to ``load``.  For
example:

.. code-block:: python

   from yt.mods import *
   filenames = ["DD0030/output_0030", "DD0040/output_0040"]
   ts = DatasetSeries.from_filenames(filenames)

This will create a new time series, populated with the output files ``DD0030``
and ``DD0040``.  This object, here called ``ts``, can now be analyzed in bulk.
Alternately, you can specify a pattern that is supplied to :mod:`glob`, and
those filenames will be sorted and returned.  Here is an example:

.. code-block:: python

   from yt.mods import *
   ts = DatasetSeries.from_filenames("*/*.index")

Analyzing Each Dataset In Sequence
----------------------------------

The :class:`~yt.data_objects.time_series.DatasetSeries` object has two primary
methods of iteration.  The first is a very simple iteration, where each object
is returned for iteration:

.. code-block:: python

   from yt.mods import *
   ts = DatasetSeries.from_filenames("*/*.index")
   for ds in ts:
       print ds.current_time

This can also operate in parallel, using
:meth:`~yt.data_objects.time_series.DatasetSeries.piter`.  For more examples,
see:

 * :ref:`parallel-time-series-analysis`
 * The cookbook recipe for :ref:`cookbook-time-series-analysis`
 * :class:`~yt.data_objects.time_series.DatasetSeries`

Prepared Time Series Analysis
-----------------------------

A few handy functions for treating time series data as a uniform, single object
are also available.

.. warning:: The future of these functions is uncertain: they may be removed in
   the future!

Simple Analysis Tasks
~~~~~~~~~~~~~~~~~~~~~

The available tasks that come built-in can be seen by looking at the output of
``ts.tasks.keys()``.  For instance, one of the simplest ones is the
``MaxValue`` task.  We can execute this task by calling it with the field whose
maximum value we want to evaluate:

.. code-block:: python

   from yt.mods import *
   ts = TimeSeries.from_filenames("*/*.index")
   max_rho = ts.tasks["MaximumValue"]("density")

When we call the task, the time series object executes the task on each
component dataset.  The results are then returned to the user.  More
complex, multi-task evaluations can be conducted by using the
:meth:`~yt.data_objects.time_series.DatasetSeries.eval` call, which accepts a
list of analysis tasks.

Analysis Tasks Applied to Objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just as some tasks can be applied to datasets as a whole, one can also apply
the creation of objects to datasets.  This means that you are able to construct
a generalized "sphere" operator that will be created inside all datasets, which
you can then calculate derived quantities (see :ref:`derived-quantities`) from.

For instance, imagine that you wanted to create a sphere that is centered on
the most dense point in the simulation and that is 1 pc in radius, and then
calculate the angular momentum vector on this sphere.  You could do that with
this script:

.. code-block:: python

   from yt.mods import *
   ts = TimeSeries.from_filenames("*/*.index")
   sphere = ts.sphere("max", (1.0, "pc"))
   L_vecs = sphere.quantities["AngularMomentumVector"]()

Note that we have specified the units differently than usual -- the time series
objects allow units as a tuple, so that in cases where units may change over
the course of several outputs they are correctly set at all times.  This script
simply sets up the time series object, creates a sphere, and then runs
quantities on it.  It is designed to look very similar to the code that would
conduct this analysis on a single output.

All of the objects listed in :ref:`available-objects` are made available in
the same manner as "sphere" was used above.

Creating Analysis Tasks
~~~~~~~~~~~~~~~~~~~~~~~

If you wanted to look at the mass in star particles as a function of time, you
would write a function that accepts params and ds and then decorate it with
analysis_task. Here we have done so:

.. code-block:: python

   @analysis_task(('particle_type',))
   def MassInParticleType(params, ds):
       dd = ds.all_data()
       ptype = (dd["particle_type"] == params.particle_type)
       return (ptype.sum(), dd["ParticleMassMsun"][ptype].sum())

   ms = ts.tasks["MassInParticleType"](4)
   print ms

This allows you to create your own analysis tasks that will be then available
to time series data objects.  Since ``DatasetSeries`` objects iterate over
filenames in parallel by default, this allows for transparent parallelization. 

.. _analyzing-an-entire-simulation:

Analyzing an Entire Simulation
------------------------------

The parameter file used to run a simulation contains all the information 
necessary to know what datasets should be available.  The ``simulation`` 
convenience function allows one to create a ``DatasetSeries`` object of all 
or a subset of all data created by a single simulation.

.. note:: Currently only implemented for Enzo.  Other simulation types coming 
   soon.

To instantiate, give the parameter file and the simulation type.

.. code-block:: python

  from yt.mods import *
  my_sim = simulation('enzo_tiny_cosmology/32Mpc_32.enzo', 'Enzo',
                      find_outputs=False)

Then, create a ``DatasetSeries`` object with the :meth:`get_time_series` 
function.  With no additional keywords, the time series will include every 
dataset.  If the **find_outputs** keyword is set to True, a search of the 
simulation directory will be performed looking for potential datasets.  These 
datasets will be temporarily loaded in order to figure out the time and 
redshift associated with them.  This can be used when simulation data was 
created in a non-standard way, making it difficult to guess the corresponding 
time and redshift information

.. code-block:: python

  my_sim.get_time_series()

After this, time series analysis can be done normally.

.. code-block:: python

  for ds in my_sim.piter()
      all_data = ds.all_data()
      print all_data.quantities['Extrema']('density')
 
Additional keywords can be given to :meth:`get_time_series` to select a subset
of the total data:

 * **time_data** (*bool*): Whether or not to include time outputs when 
   gathering datasets for time series.  Default: True.

 * **redshift_data** (*bool*): Whether or not to include redshift outputs 
   when gathering datasets for time series.  Default: True.

 * **initial_time** (*float*): The earliest time for outputs to be included.  
   If None, the initial time of the simulation is used.  This can be used in 
   combination with either final_time or final_redshift.  Default: None.

 * **final_time** (*float*): The latest time for outputs to be included.  If 
   None, the final time of the simulation is used.  This can be used in 
   combination with either initial_time or initial_redshift.  Default: None.

 * **times** (*list*): A list of times for which outputs will be found.
   Default: None.

 * **time_units** (*str*): The time units used for requesting outputs by time.
   Default: '1' (code units).

 * **initial_redshift** (*float*): The earliest redshift for outputs to be 
   included.  If None, the initial redshift of the simulation is used.  This
   can be used in combination with either final_time or final_redshift.
   Default: None.

 * **final_time** (*float*): The latest redshift for outputs to be included.  
   If None, the final redshift of the simulation is used.  This can be used 
   in combination with either initial_time or initial_redshift.  
   Default: None.

 * **redshifts** (*list*): A list of redshifts for which outputs will be found.
   Default: None.

 * **initial_cycle** (*float*): The earliest cycle for outputs to be 
   included.  If None, the initial cycle of the simulation is used.  This can
   only be used with final_cycle.  Default: None.

 * **final_cycle** (*float*): The latest cycle for outputs to be included.  
   If None, the final cycle of the simulation is used.  This can only be used 
   in combination with initial_cycle.  Default: None.

 * **tolerance** (*float*):  Used in combination with "times" or "redshifts" 
   keywords, this is the tolerance within which outputs are accepted given 
   the requested times or redshifts.  If None, the nearest output is always 
   taken.  Default: None.

 * **parallel** (*bool*/*int*): If True, the generated DatasetSeries will 
   divide the work such that a single processor works on each dataset.  If an
   integer is supplied, the work will be divided into that number of jobs.
   Default: True.

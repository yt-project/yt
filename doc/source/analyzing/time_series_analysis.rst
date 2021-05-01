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

The simplest mechanism for creating a ``DatasetSeries`` object is to pass a glob
pattern to the ``yt.load`` function.

.. code-block:: python

   import yt

   ts = yt.load("DD????/DD????")

This will create a new time series, populated with all datasets that match the
pattern "DD" followed by four digits.  This object, here called ``ts``, can now
be analyzed in bulk.  Alternately, you can specify an already formatted list of
filenames directly to the :class:`~yt.data_objects.time_series.DatasetSeries`
initializer:

.. code-block:: python

   import yt

   ts = yt.DatasetSeries(["DD0030/DD0030", "DD0040/DD0040"])

Analyzing Each Dataset In Sequence
----------------------------------

The :class:`~yt.data_objects.time_series.DatasetSeries` object has two primary
methods of iteration.  The first is a very simple iteration, where each object
is returned for iteration:

.. code-block:: python

   import yt

   ts = yt.load("*/*.index")
   for ds in ts:
       print(ds.current_time)

This can also operate in parallel, using
:meth:`~yt.data_objects.time_series.DatasetSeries.piter`.  For more examples,
see:

 * :ref:`parallel-time-series-analysis`
 * The cookbook recipe for :ref:`cookbook-time-series-analysis`
 * :class:`~yt.data_objects.time_series.DatasetSeries`

.. _analyzing-an-entire-simulation:

Analyzing an Entire Simulation
------------------------------

.. note:: Implemented for the Enzo, Gadget, OWLS, and Exodus II frontends.

The parameter file used to run a simulation contains all the information
necessary to know what datasets should be available.  The ``simulation``
convenience function allows one to create a ``DatasetSeries`` object of all
or a subset of all data created by a single simulation.

To instantiate, give the parameter file and the simulation type.

.. code-block:: python

  import yt

  my_sim = yt.load_simulation("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo")

Then, create a ``DatasetSeries`` object with the
:meth:`frontends.enzo.simulation_handling.EnzoSimulation.get_time_series`
function.  With no additional keywords, the time series will include every
dataset.  If the ``find_outputs`` keyword is set to ``True``, a search of the
simulation directory will be performed looking for potential datasets.  These
datasets will be temporarily loaded in order to figure out the time and
redshift associated with them.  This can be used when simulation data was
created in a non-standard way, making it difficult to guess the corresponding
time and redshift information

.. code-block:: python

  my_sim.get_time_series()

After this, time series analysis can be done normally.

.. code-block:: python

  for ds in my_sim.piter():
      all_data = ds.all_data()
      print(all_data.quantities.extrema("density"))

Additional keywords can be given to
:meth:`frontends.enzo.simulation_handling.EnzoSimulation.get_time_series`
to select a subset of the total data:

* ``time_data`` (*bool*): Whether or not to include time outputs when
  gathering datasets for time series.  Default: True.  (Enzo only)

* ``redshift_data`` (*bool*): Whether or not to include redshift outputs
  when gathering datasets for time series.  Default: True.  (Enzo only)

* ``initial_time`` (*float*): The earliest time for outputs to be included.
  If None, the initial time of the simulation is used.  This can be used in
  combination with either ``final_time`` or ``final_redshift``.  Default: None.

* ``final_time`` (*float*): The latest time for outputs to be included.  If
  None, the final time of the simulation is used.  This can be used in
  combination with either ``initial_time`` or ``initial_redshift``.  Default: None.

* ``times`` (*list*): A list of times for which outputs will be found.
  Default: None.

* ``initial_redshift`` (*float*): The earliest redshift for outputs to be
  included.  If None, the initial redshift of the simulation is used.  This
  can be used in combination with either ``final_time`` or ``final_redshift``.
  Default: None.

* ``final_redshift`` (*float*): The latest redshift for outputs to be included.
  If None, the final redshift of the simulation is used.  This can be used
  in combination with either ``initial_time`` or ``initial_redshift``.
  Default: None.

* ``redshifts`` (*list*): A list of redshifts for which outputs will be found.
  Default: None.

* ``initial_cycle`` (*float*): The earliest cycle for outputs to be
  included.  If None, the initial cycle of the simulation is used.  This can
  only be used with final_cycle.  Default: None.  (Enzo only)

* ``final_cycle`` (*float*): The latest cycle for outputs to be included.
  If None, the final cycle of the simulation is used.  This can only be used
  in combination with initial_cycle.  Default: None.  (Enzo only)

* ``tolerance`` (*float*):  Used in combination with ``times`` or ``redshifts``
  keywords, this is the tolerance within which outputs are accepted given
  the requested times or redshifts.  If None, the nearest output is always
  taken.  Default: None.

* ``parallel`` (*bool*/*int*): If True, the generated ``DatasetSeries`` will
  divide the work such that a single processor works on each dataset.  If an
  integer is supplied, the work will be divided into that number of jobs.
  Default: True.

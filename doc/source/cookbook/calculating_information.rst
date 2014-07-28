Calculating Dataset Information
-------------------------------

These recipes demonstrate methods of calculating quantities in a simulation,
either for later visualization or for understanding properties of fluids and
particles in the simulation.

Average Field Value
~~~~~~~~~~~~~~~~~~~

This recipe is a very simple method of calculating the global average of a
given field, as weighted by another field.

.. yt_cookbook:: average_value.py

Mass Enclosed in a Sphere
~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe constructs a sphere and then sums the total mass in particles and
fluids in the sphere.

.. yt_cookbook:: sum_mass_in_sphere.py

Global Phase Plot
~~~~~~~~~~~~~~~~~

This is a simple recipe to show how to open a dataset and then plot a couple
global phase diagrams, save them, and quit.

.. yt_cookbook:: global_phase_plots.py

Radial Velocity Profile
~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to subtract off a bulk velocity on a sphere before
calculating the radial velocity within that sphere.

.. yt_cookbook:: rad_velocity.py 

Simulation Analysis
~~~~~~~~~~~~~~~~~~~

This uses :class:`~yt.data_objects.time_series.SimulationTimeSeries` to
calculate the extrema of a series of outputs, whose names it guesses in
advance.  This will run in parallel and take advantage of multiple MPI tasks.

.. yt_cookbook:: simulation_analysis.py


.. _cookbook-time-series-analysis:

Time Series Analysis
~~~~~~~~~~~~~~~~~~~~

This recipe shows how to calculate a number of quantities on a set of parameter
files.  Note that it is parallel aware, and that if you only wanted to run in
serial the operation ``for pf in ts:`` would also have worked identically.

.. yt_cookbook:: time_series.py

.. _cookbook-simple-derived-fields:

Simple Derived Fields
~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to create a simple derived field, 
thermal_energy_density, and then generate a projection from it.

.. yt_cookbook:: derived_field.py

.. _cookbook-complex-derived-fields:

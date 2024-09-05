Calculating Dataset Information
-------------------------------

These recipes demonstrate methods of calculating quantities in a simulation,
either for later visualization or for understanding properties of fluids and
particles in the simulation.

Average Field Value
~~~~~~~~~~~~~~~~~~~

This recipe is a very simple method of calculating the global average of a
given field, as weighted by another field.
See :ref:`derived-quantities` for more information.

.. yt_cookbook:: average_value.py

Mass Enclosed in a Sphere
~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe constructs a sphere and then sums the total mass in particles and
fluids in the sphere.
See :ref:`available-objects` and :ref:`derived-quantities` for more information.

.. yt_cookbook:: sum_mass_in_sphere.py

Global Phase Plot
~~~~~~~~~~~~~~~~~

This is a simple recipe to show how to open a dataset and then plot a couple
global phase diagrams, save them, and quit.
See :ref:`how-to-make-2d-profiles` for more information.

.. yt_cookbook:: global_phase_plots.py

.. _cookbook-radial-velocity:

Radial Velocity Profile
~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to subtract off a bulk velocity on a sphere before
calculating the radial velocity within that sphere.
See :ref:`how-to-make-1d-profiles` for more information on creating profiles and
:ref:`field_parameters` for an explanation of how the bulk velocity is provided
to the radial velocity field function.

.. yt_cookbook:: rad_velocity.py

Simulation Analysis
~~~~~~~~~~~~~~~~~~~

This uses :class:`~yt.data_objects.time_series.DatasetSeries` to
calculate the extrema of a series of outputs, whose names it guesses in
advance.  This will run in parallel and take advantage of multiple MPI tasks.
See :ref:`parallel-computation` and :ref:`time-series-analysis` for more
information.

.. yt_cookbook:: simulation_analysis.py


.. _cookbook-time-series-analysis:

Time Series Analysis
~~~~~~~~~~~~~~~~~~~~

This recipe shows how to calculate a number of quantities on a set of parameter
files.  Note that it is parallel aware, and that if you only wanted to run in
serial the operation ``for pf in ts:`` would also have worked identically.
See :ref:`parallel-computation` and :ref:`time-series-analysis` for more
information.

.. yt_cookbook:: time_series.py

.. _cookbook-simple-derived-fields:

Simple Derived Fields
~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to create a simple derived field,
``thermal_energy_density``, and then generate a projection from it.
See :ref:`creating-derived-fields` and :ref:`projection-plots` for more
information.

.. yt_cookbook:: derived_field.py

.. _cookbook-complicated-derived-fields:

Complicated Derived Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to use the
:meth:`~yt.frontends.flash.data_structures.FLASHDataset.add_gradient_fields` method
to generate gradient fields and use them in a more complex derived field.

.. yt_cookbook:: hse_field.py

Using Particle Filters to Calculate Star Formation Rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe demonstrates how to use a particle filter to calculate the star
formation rate in a galaxy evolution simulation.
See :ref:`filtering-particles` for more information.

.. yt_cookbook:: particle_filter_sfr.py

Making a Turbulent Kinetic Energy Power Spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe shows how to use ``yt`` to read data and put it on a uniform
grid to interface with the NumPy FFT routines and create a turbulent
kinetic energy power spectrum.  (Note: the dataset used here is of low
resolution, so the turbulence is not very well-developed.  The spike
at high wavenumbers is due to non-periodicity in the z-direction).

.. yt_cookbook:: power_spectrum_example.py

Downsampling an AMR Dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This recipe shows how to use the ``max_level`` attribute of a yt data
object to only select data up to a maximum AMR level.

.. yt_cookbook:: downsampling_amr.py

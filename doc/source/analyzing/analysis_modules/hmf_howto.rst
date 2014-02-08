.. _hmf_howto:

Halo Mass Function: Start to Finish
===================================

This how-to steps through the three simple steps it takes to find the halo
mass function for an ``enzo`` dataset. The first step is to find the haloes,
second to find the virial mass (gas + dark matter) contained in each halo using
the halo profiler, and
finally build a halo mass function of the haloes and an analytical fit
for comparison.

Halo Finding
------------

The first step is find the haloes in the simulation. There are a number of ways
to do this explained in detail in :ref:`halo_finding`.
You are encouraged to read about the differences between the halo finders and
experiment with different settings and functions.
For the purposes of
the halo mass function, Friends of Friends (FOF) is probably not the best choice.
Therefore, below
is a simple example of how to run HOP on a dataset. This example will also
write out a text file with the halo particulars, which will be used in the next
step.

.. code-block:: python

  from yt.mods import *
  pf = load("data0001")
  halo_list = HaloFinder(pf)
  halo_list.write_out("HopAnalysis.out")

The only important columns of data in the text file ``HopAnalysis.out``
are the halo number, center of
mass, and maximum radius. These are the only values that the next step requires.

Halo Profiling
--------------

The halo profiler (:ref:`halo_profiling`) is a powerful tool that can analyze
haloes in many ways. It is beneficial to read its documentation to become
familiar with it before using it.
For this exercise, only the virial mass of each
halo is important. This script below will take the output from the previous step
and find the virial mass for each halo, and save it to a text file.

.. code-block:: python

  import yt.analysis_modules.halo_profiler.api as HP
  hp = HP.HaloProfiler("data0001", halo_list_file='HopAnalysis.out')
  hp.add_halo_filter(HP.VirialFilter,must_be_virialized=True,
                overdensity_field='ActualOverdensity',
                virial_overdensity=200,
                virial_filters=[['TotalMassMsun','>=','1e8']],
                virial_quantities=['TotalMassMsun','RadiusMpc'])
  hp.make_profiles(filename="VirialHaloes.out")

This script limits the output to virialized haloes with mass greater than or
equal to 1e8 solar masses. If you run into problems, try pre-filtering problem
haloes (:ref:`halo_profiler_pre_filters`).

Halo Mass Function
------------------

The halo mass function extension (see :ref:`halo_mass_function`) reads in the
contents of ``VirialHaloes.out`` and can output files which can be used
to make a plot. In this script below, the results from the previous step are
read in, analyzed, and at the same time an analytical fit is calculated, and
finally the results are written out to two files. The plot found from the haloes
is saved to a file named ``hmf-haloes.dat`` and the fit to ``hmf-fit.dat``.
For these files, the columns that are most likely to be needed for plotting are
the first (halo mass bin), and the fourth and third (cumulative number density
of haloes) for the ``-fit.dat`` and ``-haloes.dat`` files. See
the full halo mass function documentation for more details on the contents of
the files.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_mass_function.api import *
  pf = load("data0001")
  hmf = HaloMassFcn(pf, halo_file="VirialHaloes.out", 
  sigma8input=0.9, primordial_index=1., omega_baryon0=0.06,
  fitting_function=4, mass_column=5, num_sigma_bins=200)
  hmf.write_out(prefix='hmf')

Inside the call to ``HaloMassFcn`` there are several hand-set parameters that 
*must* be correct for the analytical fit to be correct. The three cosmological
parameters (``sigma8input``, ``primordial_index`` and ``omega_baryon0``) are
not stored with ``Enzo`` datasets, so they must be found from the initial
conditions of the simulation. ``mass_column`` is set to 5 because that is the
zero-indexed ordinal of the mass column in the ``VirialHaloes.out`` file.
``num_sigma_bins`` controls how many mass bins the haloes are dropped into,
and ``fitting_function`` controls which analytical function to use.

Putting it All Together
-----------------------

It is not necessary to run each step separately from the others. This example
below will run all steps at once.

.. code-block:: python

  from yt.mods import *
  import yt.analysis_modules.halo_profiler.api as HP
  from yt.analysis_modules.halo_mass_function.api import *
  
  # If desired, start loop here.
  pf = load("data0001")
  
  halo_list = HaloFinder(pf)
  halo_list.write_out("HopAnalysis.out")
  
  hp = HP.HaloProfiler("data0001", halo_list_file='HopAnalysis.out')
  hp.add_halo_filter(HP.VirialFilter,must_be_virialized=True,
                overdensity_field='ActualOverdensity',
                virial_overdensity=200,
                virial_filters=[['TotalMassMsun','>=','1e8']],
                virial_quantities=['TotalMassMsun','RadiusMpc'])
  hp.make_profiles(filename="VirialHaloes.out")
  
  hmf = HaloMassFcn(pf, halo_file="VirialHaloes.out", 
  sigma8input=0.9, primordial_index=1., omega_baryon0=0.06,
  fitting_function=4, mass_column=5, num_sigma_bins=200)
  hmf.write_out(prefix='hmf')
  # End loop here.

The script above will work in parallel which can reduce runtimes substantially.
If this analysis is to be run on a sequence of datasets, the section that needs
to be inside the loop is shown bracketed by comments. Be careful how the
output files are named as to not over-write output from previous loop cycles.

Plotting
--------

When plotting the output, be careful about the units of the output for the
halo mass function. The figure shown in the documentation (on this page:
:ref:`halo_mass_function`) has the number density of haloes per (h^-1 Mpc)^3,
which is different than the output of the halo mass extension (which is
haloes per (Mpc)^3). To get the same units as the figure for the ``-fit.dat``
and ``-haloes.dat`` files, divide the fourth and third column by the comoving
volume cubed of the simulation, respectively when plotting.

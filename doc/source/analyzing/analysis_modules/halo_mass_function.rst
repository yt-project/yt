.. _halo_mass_function:

Halo Mass Function
==================

The Halo Mass Function extension is capable of outputting the halo mass function
for a collection halos (input), and/or an analytical fit over a given mass range
for a set of specified cosmological parameters.
This extension is based on code generously provided by Brian O'Shea.

General Overview
----------------

A halo mass function can be created for the halos identified in a cosmological 
simulation, as well as analytic fits using any arbitrary set of cosmological
paramters. In order to create a mass function for simulated halos, they must
first be identified (using HOP, FOF, or Rockstar, see 
:ref:`halo_catalog`) and loaded as a halo dataset object. The distribution of
halo masses will then be found, and can be compared to the analytic prediction
at the same redshift and using the same cosmological parameters as were used
in the simulation. Care should be taken in this regard, as the analytic fit
requires the specification of cosmological parameters that are not necessarily 
stored in the halo or simulation datasets, and must be specified by the user.
Efforts have been made to set reasonable defaults for these parameters, but 
setting them to identically match those used in the simulation will produce a
much better comparison.

Analytic halo mass functions can also be created without a halo dataset by 
providing either a simulation dataset or specifying cosmological parameters by
hand. yt includes 5 analytic fits for the halo mass function which can be
selected.


Analytical Fits
---------------

There are five analytical fits to choose from.

  1. `Press-Schechter (1974) <http://adsabs.harvard.edu/abs/1974ApJ...187..425P>`_
  2. `Jenkins (2001) <http://adsabs.harvard.edu/abs/2001MNRAS.321..372J>`_
  3. `Sheth-Tormen (2002) <http://adsabs.harvard.edu/abs/2002MNRAS.329...61S>`_
  4. `Warren (2006) <http://adsabs.harvard.edu/abs/2006ApJ...646..881W>`_
  5. `Tinker (2008) <http://adsabs.harvard.edu/abs/2008ApJ...688..709T>`_

We encourage reading each of the primary sources.
In general, we recommend the Warren fitting function because it matches
simulations over a wide range of masses very well.
The Warren fitting function is the default (equivalent to not specifying
``fitting_function`` in ``HaloMassFcn()``, below).
The Tinker fit is for the :math:`\Delta=300` fits given in the paper, which
appears to fit HOP threshold=80.0 fairly well.


Basic Halo Mass Function Creation
---------------------------------

The simplest way to create a halo mass function object is to simply pass it no
arguments and let it use the default cosmological parameters.

.. code-block:: python

  from yt.analysis_modules.halo_mass_function.api import *

  hmf = HaloMassFcn()

This will create a HaloMassFcn object off of which arrays holding the information
about the analytic mass function hang. Creating the halo mass function for a set
of simulated halos requires only the loaded halo dataset to be passed as an 
argument. This also creates the analytic mass function using all parameters that 
can be extracted from the halo dataset, at the same redshift, spanning a similar
range of halo masses.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_mass_function.api import *

  my_halos = load("rockstar_halos/halos_0.0.bin")
  hmf = HaloMassFcn(halos_ds=my_halos)

A simulation dataset can be passed along with additonal cosmological parameters 
to create an analytic mass function.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_mass_function.api import *

  my_ds = load("RD0027/RedshiftOutput0027")
  hmf = HaloMassFcn(simulation_ds=my_ds, omega_baryon0=0.05, primordial_index=0.96, 
                    sigma8 = 0.8, log_mass_min=5, log_mass_max=9)

The analytic mass function can be created for a set of arbitrary cosmological 
parameters without any dataset being passed as an argument.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_mass_function.api import *

  hmf = HaloMassFcn(omega_baryon0=0.05, omega_matter0=0.27, 
                    omega_lambda0=0.73, hubble0=0.7, this_redshift=10,
                    log_mass_min=5, log_mass_max=9, fitting_function=5)

Keyword Arguments
-----------------

* **simulation_ds** (*Simulation dataset object*)
  The loaded simulation dataset, used to set cosmological paramters.
  Default : None.

* **halos_ds** (*Halo dataset object*)
  The halos from a simulation to be used for creation of the 
  halo mass function in the simulation.
  Default : None.

* **make_analytic** (*bool*)
  Whether or not to calculate the analytic mass function to go with 
  the simulated halo mass function.  Automatically set to true if a 
  simulation dataset is provided.
  Default : True.

* **omega_matter0** (*float*)
  The fraction of the universe made up of matter (dark and baryonic). 
  Default : 0.2726.

* **omega_lambda0** (*float*)
  The fraction of the universe made up of dark energy. 
  Default : 0.7274.

* **omega_baryon0**  (*float*)
  The fraction of the universe made up of baryonic matter. This is not 
  always stored in the datset and should be checked by hand.
  Default : 0.0456.

* **hubble0** (*float*)
  The expansion rate of the universe in units of 100 km/s/Mpc. 
  Default : 0.704.

* **sigma8** (*float*)
  The amplitude of the linear power spectrum at z=0 as specified by 
  the rms amplitude of mass-fluctuations in a top-hat sphere of radius 
  8 Mpc/h. This is not always stored in the datset and should be 
  checked by hand.
  Default : 0.86.

* **primoridal_index** (*float*)
  This is the index of the mass power spectrum before modification by 
  the transfer function. A value of 1 corresponds to the scale-free 
  primordial spectrum. This is not always stored in the datset and 
  should be checked by hand.
  Default : 1.0.

* **this_redshift** (*float*)
  The current redshift. 
  Default : 0.

* **log_mass_min** (*float*)
  The log10 of the mass of the minimum of the halo mass range. This is
  set automatically by the range of halo masses if a simulated halo 
  dataset is provided. If a halo dataset if not provided and no value
  is specified, it will be set to 5. Units: M_solar
  Default : None.

* **log_mass_max** (*float*)
  The log10 of the mass of the maximum of the halo mass range. This is
  set automatically by the range of halo masses if a simulated halo 
  dataset is provided. If a halo dataset if not provided and no value
  is specified, it will be set to 16. Units: M_solar
  Default : None.

* **num_sigma_bins** (*float*)
  The number of bins (points) to use for the calculation of the 
  analytic mass function. 
  Default : 360.

* **fitting_function** (*int*)
  Which fitting function to use. 1 = Press-Schechter, 2 = Jenkins, 
  3 = Sheth-Tormen, 4 = Warren, 5 = Tinker
  Default : 4.

Outputs
-------

A HaloMassFnc object has several arrays hanging off of it containing the 

* **masses_sim**: Halo masses from simulated halos. Units: M_solar

* **n_cumulative_sim**: Number density of halos with mass greater than the 
  corresponding mass in masses_sim. Units: comoving Mpc^-3

* **masses_analytic**: Masses used for the generation of the analytic mass 
  function. Units: M_solar

* **n_cumulative_analytic**: Number density of halos with mass greater then 
  the corresponding mass in masses_analytic. Units: comoving Mpc^-3

* **dndM_dM_analytic**: Differential number density of halos, (dn/dM)*dM.

After the mass function has been created for both simulated halos and the
corresponding analytic fits, they can be plotted though something along the 
lines of

.. code-block:: python

  import yt
  from yt.analysis_modules.halo_mass_function.api import *
  import matplotlib.pyplot as plt

  my_halos = yt.load("rockstar_halos/halos_0.0.bin")
  hmf = HaloMassFcn(halos_ds=my_halos)

  plt.loglog(hmf.masses_sim, hmf.n_cumulative_sim)
  plt.loglog(hmf.masses_analytic, hmf.n_cumulative_analytic)

Attached to ``hmf`` is the convenience function ``write_out``, which saves the 
halo mass function to a text file. (continued from above)
.. code-block:: python

  hmf.write_out(prefix='hmf', analytic=True, simulated=True)

This writes the files ``hmf-analytic.dat`` with columns:

* mass [Msun]
* cumulative number density of halos [comoving Mpc^-3]
* (dn/dM)*dM (differential number density of halos) [comoving Mpc^-3]

and the file ``hmf-simulated.dat`` with columns:

* mass [Msun]
* cumulative number density of halos [comoving Mpc^-3]

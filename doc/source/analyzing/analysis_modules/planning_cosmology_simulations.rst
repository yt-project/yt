.. _planning-cosmology-simulations:

Planning Simulations to use LightCones or LightRays
===================================================

If you want to run a cosmological simulation that will have just enough data 
outputs to create a cosmology splice, the 
:meth:`~yt.analysis_modules.cosmological_observation.cosmology_splice.CosmologySplice.plan_cosmology_splice` 
function will calculate a list of redshifts outputs that will minimally 
connect a redshift interval.

.. code-block:: python

  from yt.analysis_modules.cosmological_observation.api import CosmologySplice
  my_splice = CosmologySplice('enzo_tiny_cosmology/32Mpc_32.enzo', 'Enzo')
  my_splice.plan_cosmology_splice(0.0, 0.1, filename='redshifts.out')

This will write out a file, formatted for simulation type, with a list of 
redshift dumps.  The keyword arguments are:

* ``decimals`` (*int*): The decimal place to which the output redshift will 
  be rounded.  If the decimal place in question is nonzero, the redshift will 
  be rounded up to ensure continuity of the splice.  Default: 3.

* ``filename`` (*str*): If provided, a file will be written with the redshift 
  outputs in the form in which they should be given in the enzo parameter 
  file.  Default: None.

* ``start_index`` (*int*): The index of the first redshift output.  Default: 0.

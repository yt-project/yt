.. _comoving_units_and_code_units:

Comoving units and code units
=============================

.. notebook:: 3)_Comoving_units_and_code_units.ipynb

.. _cosmological-units:

Units for Cosmological Datasets
-------------------------------

yt has additional capabilities to handle the comoving coordinate system used
internally in cosmological simulations. Simulations that use comoving
coordinates, all length units have three other counterparts correspoding to
comoving units, scaled comoving units, and scaled proper units. In all cases
'scaled' units refer to scaling by the reduced Hubble parameter - i.e. the length
unit is what it would be in a universe where Hubble's parameter is 100 km/s/Mpc.

To access these different units, yt has a common naming system. Scaled units are denoted by
dividing by the scaled Hubble parameter ``h`` (which is in itself a unit). Comoving
units are denoted by appending ``cm`` to the end of the unit name.

Using the parsec as an example,

``pc``
    Proper parsecs, :math:`\rm{pc}`.

``pccm``
    Comoving parsecs, :math:`\rm{pc}/(1+z)`.

``pccm/h``
    Comoving parsecs normalized by the scaled hubble constant, :math:`\rm{pc}/h/(1+z)`.

``pc/h``
    Proper parsecs, normalized by the scaled hubble constant, :math:`\rm{pc}/h`.

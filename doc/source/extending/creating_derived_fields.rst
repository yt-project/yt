Creating Derived Fields
=======================

One of the more powerful means of extending yt is through the usage of derived
fields.  These are fields that describe a value at each cell in a simulation.

Defining a New Field
--------------------

So once a new field has been 

Field Options
-------------

How Do Units Work?
------------------

Everything is done under the assumption that all of the native Enzo fields that
yt knows about are converted to cgs before being handed to any processing
routines.

Which Enzo Fields Does yt Know About?
-------------------------------------

* Density
* Temperature
* Gas Energy
* Total Energy
* [xyz]-velocity
* Species fields: HI, HII, Electron, HeI, HeII, HeIII, HM, H2I, H2II, DI, DII, HDI
* Particle mass, velocity, 
=======

Derived fields are one of the more powerful ways in which you can extend yt.

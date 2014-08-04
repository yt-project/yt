.. _ellipsoid_analysis:

Halo Ellipsoid Analysis
=======================
.. sectionauthor:: Geoffrey So <gso@physics.ucsd.edu>

.. warning:: This is my first attempt at modifying the yt source code,
   so the program may be bug ridden.  Please send yt-dev an email and
   address to Geoffrey So if you discover something wrong with this
   portion of the code.

Purpose
-------

The purpose of creating this feature in yt is to analyze field
properties that surround dark matter haloes.  Originally, this was
usually done with the sphere 3D container, but since many halo
particles are linked together in a more elongated shape, I thought it
would be better to use an ellipsoid 3D container to wrap around the
particles.  This way, less of the empty-of-particle space around the
halo would be included when doing the analysis of field properties
where the particles are suppose to occupy.

General Overview
----------------

In order to use the ellipsoid 3D container object, one must supply it
with a center, the magnitude of the semi-principle axes, the direction
of the first semi-principle axis, the tilt angle (rotation angle about
the y axis that will align the first semi-principle axis with the x
axis once it is aligned in the x-z plane.)

Once those parameters are determined, the function "ellipsoid" will
return the 3D object, and users will be able to get field attributes
from the data object just as they would from spheres, cylinders etc.

Example
-------

To use the ellipsoid container to get field information, you
will have to first determine the ellipsoid's parameters.  This can be
done with the haloes obtained from halo finding, but essentially it
takes the information:

  #. Center position x,y,z
  #. List of particles position x,y,z

And calculates the ellipsoid information needed for the 3D container.

What I usually do is get this information from the halo finder output
files in the .h5 HDF5 binary format. I load them into memory using the
LoadHaloes() function instead of reading in the ASCII output.

Halo Finding
~~~~~~~~~~~~
.. code-block:: python

  import yt
  from yt.analysis_modules.halo_finding.api import *

  ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
  halo_list = parallelHF(ds)
  halo_list.dump('MyHaloList')

Ellipsoid Parameters
~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

  import yt
  from yt.analysis_modules.halo_finding.api import *

  ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
  haloes = LoadHaloes(ds, 'MyHaloList')

Once the halo information is saved you can load it into the data
object "haloes", you can get loop over the list of haloes and do

.. code-block:: python

  ell_param = haloes[0].get_ellipsoid_parameters()

This will return 6 items

#. The center of mass as an array.
#. A as a float.  (Must have A>=B)
#. B as a float.  (Must have B>=C)
#. C as a float.  (Must have C > cell size)
#. e0 vector as an array.  (now normalized automatically in the code)
#. tilt as a float.

The center of mass would be the same one as returned by the halo
finder.  The A, B, C are the largest to smallest magnitude of the
ellipsoid's semi-principle axes. "e0" is the largest semi-principle
axis vector direction that would have magnitude A but normalized.  
The "tilt" is an angle measured in radians.  It can be best described
as after the rotation about the z-axis to allign e0 to x in the x-y
plane, and then rotating about the y-axis to align e0 completely to
the x-axis, the angle remaining to rotate about the x-axis to align
both e1 to the y-axis and e2 to the z-axis.

Ellipsoid 3D Container
~~~~~~~~~~~~~~~~~~~~~~

Once the parameters are obtained from the get_ellipsoid_parameters()
function, or picked at random by the user, it can be input into the
ellipsoid container as:

.. code-block:: python

  ell = ds.ellipsoid(ell_param[0],
  ell_param[1],
  ell_param[2],
  ell_param[3],
  ell_param[4],
  ell_param[5])
  dens = ell.quantities['TotalQuantity']('density')[0]

This way, "ell" will be the ellipsoid container, and "dens" will be
the total density of the ellipsoid in an unigrid simulation.  One can
of course use this container object with parameters that they come up
with, the ellipsoid parameters do not have to come from the Halo
Finder.  And of course, one can use the ellipsoid container with other
derived fields or fields that they are interested in.

Drawbacks
---------

Since this is a first attempt, there are many drawbacks and corners
cut.  Many things listed here will be amended when I have time.

* The ellipsoid 3D container like the boolean object, do not contain 
  particle position and velocity information.
* This currently assume periodic boundary condition, so if an
  ellipsoid center is at the edge, it will return part of the opposite
  edge field information.  Will try to put in the option to turn off
  periodicity in the future.
* This method gives a minimalistic ellipsoid centered around the
  center of mass that contains all the particles, but sometimes people
  prefer an inertial tensor triaxial ellipsoid described in 
  `Dubinski, Carlberg 1991
  <http://adsabs.harvard.edu/abs/1991ApJ...378..496D>`_.  I have that
  method composed but it is not fully tested yet.
* The method to obtain information from the halo still uses the center
  of mass as the center of the ellipsoid, so it is not making the
  smallest ellipsoid that contains the particles as possible.  To
  start at the center of the particles based on position will require
  an O(:math:`N^2`) operation, right now I'm trying to limit
  everything to O(:math:`N`) operations.  If particle count does not
  get too large, I may implement the O(:math:`N^2`) operation.
* Currently the list of haloes can be analyzed using object
  parallelism (one halo per core), but I'm not sure if haloes will get
  big enough soon that other forms of parallelism will be needed to
  analyze them due to memory constraint.
* This has only been tested on unigrid simulation data, not AMR.  In
  unigrid simulations, I can take "dens" from the example and divide
  it by the total number of cells to get the average density, in AMR
  one would need to do an volume weighted average instead.

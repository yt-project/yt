.. _astropy-integrations:

AstroPy Integrations
====================

yt enables a number of integrations with the AstroPy package. These
are listed below, but more detailed descriptions are given at the
given documentation links.

Round-Trip Unit Conversions Between yt and AstroPy
--------------------------------------------------


FITS Image File Writing
-----------------------

Converting Field Container and 1D Profile Data to AstroPy Tables
----------------------------------------------------------------

Data in field containers, such as spheres, rectangular regions, rays, 
cylinders, etc., are represented as 1D YTArrays. A set of these arrays
can then be exported to an 
`AstroPy Table <http://docs.astropy.org/en/stable/table/>`_ object, 
specifically a 
`QTable <http://docs.astropy.org/en/stable/table/mixin_columns.html#quantity-and-qtable>`_.
``QTable``s are unit-aware, and can be manipulated in a number of ways
and written to disk in several formats, including ASCII text or FITS 
files. For more details, see :ref:`fields-astropy-export`. 

Similarly, 1D profile objects can also be exported to AstroPy 
``QTable``s, optionally writing all of the profile bins or only the ones
which are used. For more details, see :ref:`profile-astropy-export`.

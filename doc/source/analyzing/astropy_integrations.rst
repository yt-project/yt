.. _astropy-integrations:

AstroPy Integrations
====================

yt enables a number of integrations with the AstroPy package. These
are listed below, but more detailed descriptions are given at the
given documentation links.

Round-Trip Unit Conversions Between yt and AstroPy
--------------------------------------------------

AstroPy has a `symbolic units implementation <https://docs.astropy.org/en/stable/units/>`_
similar to that in yt. For this reason, we have implemented "round-trip"
conversions between :class:`~yt.units.yt_array.YTArray` objects
and AstroPy's :class:`~astropy.units.Quantity` objects. These are implemented
in the :meth:`~yt.units.yt_array.YTArray.from_astropy` and
:meth:`~yt.units.yt_array.YTArray.to_astropy` methods.

FITS Image File Reading and Writing
-----------------------------------

Reading and writing FITS files is supported in yt using
`AstroPy's FITS file handling. <https://docs.astropy.org/en/stable/io/fits/>`_

yt has basic support for reading two and three-dimensional image data from FITS
files. Some limited ability to parse certain types of data (e.g., spectral cubes,
images with sky coordinates, images written using the
:class:`~yt.visualization.fits_image.FITSImageData` class described below) is
possible. See :ref:`loading-fits-data` for more information.

Fixed-resolution two-dimensional images generated from datasets using yt (such as
slices or projections) and fixed-resolution three-dimensional grids can be written
to FITS files using yt's :class:`~yt.visualization.fits_image.FITSImageData` class
and its subclasses. Multiple images can be combined into a single file, operations
can be performed on the images and their coordinates, etc. See :ref:`writing_fits_images`
for more information.

Converting Field Container and 1D Profile Data to AstroPy Tables
----------------------------------------------------------------

Data in field containers, such as spheres, rectangular regions, rays,
cylinders, etc., are represented as 1D YTArrays. A set of these arrays
can then be exported to an
`AstroPy Table <http://docs.astropy.org/en/stable/table/>`_ object,
specifically a
`QTable <http://docs.astropy.org/en/stable/table/mixin_columns.html#quantity-and-qtable>`_.
``QTable`` is unit-aware, and can be manipulated in a number of ways
and written to disk in several formats, including ASCII text or FITS
files.

Similarly, 1D profile objects can also be exported to AstroPy
``QTable``, optionally writing all of the profile bins or only the ones
which are used. For more details, see :ref:`profile-astropy-export`.

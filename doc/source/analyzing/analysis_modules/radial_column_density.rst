.. _radial-column-density:

Radial Column Density
=====================
.. sectionauthor:: Stephen Skory <s@skory.us>
.. versionadded:: 2.3

.. note:: 

    As of :code:`yt-3.0`, the radial column density analysis module is not
    currently functional.  This functionality is still available in
    :code:`yt-2.x`.  If you would like to use these features in :code:`yt-3.x`,
    help is needed to port them over.  Contact the yt-users mailing list if you
    are interested in doing this.

This module allows the calculation of column densities around a point over a
field such as ``NumberDensity`` or ``Density``.
This uses :ref:`healpix_volume_rendering` to interpolate column densities
on the grid cells.

Details
-------

This module allows the calculation of column densities around a single point.
For example, this is useful for looking at the gas around a radiating source.
Briefly summarized, the calculation is performed by first creating a number
of HEALPix shells around the central point.
Next, the value of the column density at cell centers is found by
linearly interpolating the values on the inner and outer shell.
This is added as derived field, which can be used like any other derived field.

Basic Example
-------------

In this simple example below, the radial column density for the field
``NumberDensity`` is calculated and added as a derived field named
``RCDNumberDensity``.
The calculations will use the starting point of (x, y, z) = (0.5, 0.5, 0.5) and
go out to a maximum radius of 0.5 in code units.
Due to the way normalization is handled in HEALPix, the column density
calculation can extend out only as far as the nearest face of the volume.
For example, with a center point of (0.2, 0.3, 0.4), the column density
is calculated out to only a radius of 0.2.
The column density will be output as zero (0.0) outside the maximum radius.
Just like a real number column density, when the derived is added using
``add_field``, we give the units as :math:`1/\rm{cm}^2`.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.radial_column_density.api import *
  ds = load("data0030")
  
  rcdnumdens = RadialColumnDensity(ds, 'NumberDensity', [0.5, 0.5, 0.5],
    max_radius = 0.5)
  def _RCDNumberDensity(field, data, rcd = rcdnumdens):
      return rcd._build_derived_field(data)
  add_field('RCDNumberDensity', _RCDNumberDensity, units=r'1/\rm{cm}^2')
  
  dd = ds.all_data()
  print dd['RCDNumberDensity']

The field ``RCDNumberDensity`` can be used just like any other derived field
in yt.

Additional Parameters
---------------------

Each of these parameters is added to the call to ``RadialColumnDensity()``,
just like ``max_radius`` is used above.

  * ``steps`` : integer - Because this implementation uses linear
    interpolation to calculate the column
    density at each cell, the accuracy of the solution goes up as the number of
    HEALPix surfaces is increased.
    The ``steps`` parameter controls the number of HEALPix surfaces, and a larger
    number is more accurate, but slower. Default = 10.

  * ``base`` : string - This controls where the surfaces are placed, with
    linear "lin" or logarithmic "log" spacing. The inner-most
    surface is always set to the size of the smallest cell.
    Default = "lin". 

  * ``Nside`` : int
    The resolution of column density calculation as performed by
    HEALPix. Higher numbers mean higher quality. Max = 8192.
    Default = 32.

  * ``ang_divs`` : imaginary integer
    This number controls the gridding of the HEALPix projection onto
    the spherical surfaces. Higher numbers mean higher quality.
    Default = 800j.


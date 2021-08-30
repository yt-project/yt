.. _generating-processed-data:

Generating Processed Data
=========================

Although yt provides a number of built-in visualization methods that can
process data and construct from that plots, it is often useful to generate the
data by hand and construct plots which can then be combined with other plots,
modified in some way, or even (gasp) created and modified in some other tool or
program.

.. _exporting-container-data:

Exporting Container Data
------------------------

Fields from data containers such as regions, spheres, cylinders, etc. can be exported
tabular format using either a :class:`~pandas.DataFrame` or an :class:`~astropy.table.QTable`.

To export to a :class:`~pandas.DataFrame`, use
:meth:`~yt.data_objects.data_containers.YTDataContainer.to_dataframe`:

.. code-block:: python

    sp = ds.sphere("c", (0.2, "unitary"))
    df2 = sp.to_dataframe([("gas", "density"), ("gas", "temperature")])

To export to a :class:`~astropy.table.QTable`, use
:meth:`~yt.data_objects.data_containers.YTDataContainer.to_astropy_table`:

.. code-block:: python

    sp = ds.sphere("c", (0.2, "unitary"))
    at2 = sp.to_astropy_table(fields=[("gas", "density"), ("gas", "temperature")])

For exports to :class:`~pandas.DataFrame` objects, the unit information is lost, but for
exports to :class:`~astropy.table.QTable` objects, the :class:`~yt.units.yt_array.YTArray`
objects are converted to :class:`~astropy.units.Quantity` objects.

.. _generating-2d-image-arrays:

2D Image Arrays
---------------

When making a slice, a projection or an oblique slice in yt, the resultant
:class:`~yt.data_objects.data_containers.YTSelectionContainer2D` object is created and
contains flattened arrays of the finest available data.  This means a set of
arrays for the x, y, (possibly z), dx, dy, (possibly dz) and data values, for
every point that constitutes the object.


This presents something of a challenge for visualization, as it will require
the transformation of a variable mesh of points consisting of positions and
sizes into a fixed-size array that appears like an image.  This process is that
of pixelization, which yt handles transparently internally.  You can access
this functionality by constructing a
:class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer` and supplying
to it your :class:`~yt.data_objects.data_containers.YTSelectionContainer2D`
object, as well as some information about how you want the final image to look.
You can specify both the bounds of the image (in the appropriate x-y plane) and
the resolution of the output image.  You can then have yt pixelize any field
you like.

.. note:: In previous versions of yt, there was a special class of
          FixedResolutionBuffer for off-axis slices.  This is no longer
          necessary.

To create :class:`~yt.data_objects.data_containers.YTSelectionContainer2D` objects, you can
access them as described in :ref:`data-objects`, specifically the section
:ref:`available-objects`.  Here is an example of how to window into a slice
of resolution(512, 512) with bounds of (0.3, 0.5) and (0.6, 0.8).  The next
step is to generate the actual 2D image array, which is accomplished by
accessing the desired field.

.. code-block:: python

   sl = ds.slice(0, 0.5)
   frb = FixedResolutionBuffer(sl, (0.3, 0.5, 0.6, 0.8), (512, 512))
   my_image = frb["density"]

This image may then be used in a hand-constructed Matplotlib image, for instance using
:func:`~matplotlib.pyplot.imshow`.

The buffer arrays can be saved out to disk in either HDF5 or FITS format:

.. code-block:: python

   frb.save_as_dataset("my_images.h5", fields=[("gas", "density"), ("gas", "temperature")])
   frb.export_fits(
       "my_images.fits",
       fields=[("gas", "density"), ("gas", "temperature")],
       clobber=True,
       units="kpc",
   )

In the HDF5 case, the created file can be reloaded just like a regular dataset with
``yt.load`` and will, itself, be a first-class dataset.  For more information on
this, see :ref:`saving-grid-data-containers`.
In the FITS case, there is an option for setting the ``units`` of the coordinate system in
the file. If you want to overwrite a file with the same name, set ``clobber=True``.

The :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer` can even be exported
as a 2D dataset itself, which may be operated on in the same way as any other dataset in yt:

.. code-block:: python

   ds_frb = frb.export_dataset(
       fields=[("gas", "density"), ("gas", "temperature")], nprocs=8
   )
   sp = ds_frb.sphere("c", (100.0, "kpc"))

where the ``nprocs`` parameter can be used to decompose the image into ``nprocs`` number of grids.

.. _generating-profiles-and-histograms:

Profiles and Histograms
-----------------------

Profiles and histograms can also be generated using the
:class:`~yt.visualization.profile_plotter.ProfilePlot` and
:class:`~yt.visualization.profile_plotter.PhasePlot` functions
(described in :ref:`how-to-make-1d-profiles` and
:ref:`how-to-make-2d-profiles`).  These generate profiles transparently, but the
objects they handle and create can be handled manually, as well, for more
control and access.  The :func:`~yt.data_objects.profiles.create_profile` function
can be used to generate 1, 2, and 3D profiles.

Profile objects can be created from any data object (see :ref:`data-objects`,
specifically the section :ref:`available-objects` for more information) and are
best thought of as distribution calculations.  They can either sum up or average
one quantity with respect to one or more other quantities, and they do this over
all the data contained in their source object.  When calculating average values,
the standard deviation will also be calculated.

To generate a profile, one need only specify the binning fields and the field
to be profiled.  The binning fields are given together in a list.  The
:func:`~yt.data_objects.profiles.create_profile` function will guess the
dimensionality of the profile based on the number of fields given.  For example,
a one-dimensional profile of the mass-weighted average temperature as a function of
density within a sphere can be created in the following way:

.. code-block:: python

   import yt

   ds = yt.load("galaxy0030/galaxy0030")
   source = ds.sphere("c", (10, "kpc"))
   profile = source.profile(
       [("gas", "density")],  # the bin field
       [
           ("gas", "temperature"),  # profile field
           ("gas", "radial_velocity"),
       ],  # profile field
       weight_field=("gas", "mass"),
   )

The binning, weight, and profile data can now be access as:

.. code-block:: python

   print(profile.x)  # bin field
   print(profile.weight)  # weight field
   print(profile["gas", "temperature"])  # profile field
   print(profile["gas", "radial_velocity"])  # profile field

The ``profile.used`` attribute gives a boolean array of the bins which actually
have data.

.. code-block:: python

   print(profile.used)

If a weight field was given, the profile data will represent the weighted mean
of a field.  In this case, the weighted standard deviation will be calculated
automatically and can be access via the ``profile.standard_deviation``
attribute.

.. code-block:: python

   print(profile.standard_deviation["gas", "temperature"])

A two-dimensional profile of the total gas mass in bins of density and
temperature can be created as follows:

.. code-block:: python

   profile2d = source.profile(
       [
           ("gas", "density"),
           ("gas", "temperature"),
       ],  # the x bin field  # the y bin field
       [("gas", "mass")],  # the profile field
       weight_field=None,
   )

Accessing the x, y, and profile fields work just as with one-dimensional profiles:

.. code-block:: python

   print(profile2d.x)
   print(profile2d.y)
   print(profile2d["gas", "mass"])

One of the more interesting things that is enabled with this approach is
the generation of 1D profiles that correspond to 2D profiles.  For instance, a
phase plot that shows the distribution of mass in the density-temperature
plane, with the average temperature overplotted.  The
:func:`~matplotlib.pyplot.pcolormesh` function can be used to manually plot
the 2D profile.  If you want to generate a default profile plot, you can simply
call:::

  profile.plot()

Three-dimensional profiles can be generated and accessed following
the same procedures.  Additional keyword arguments are available to control
the following for each of the bin fields: the number of bins, min and max, units,
whether to use a log or linear scale, and whether or not to do accumulation to
create a cumulative distribution function.  For more information, see the API
documentation on the :func:`~yt.data_objects.profiles.create_profile` function.

For custom bins the other keyword arguments can be overridden using the
``override_bins`` keyword argument. This accepts a dictionary with an array
for each bin field or ``None`` to use the default settings.

.. code-block:: python

    custom_bins = np.array([1e-27, 1e-25, 2e-25, 5e-25, 1e-23])
    profile2d = source.profile(
        [("gas", "density"), ("gas", "temperature")],
        [("gas", "mass")],
        override_bins={("gas", "density"): custom_bins, ("gas", "temperature"): None},
    )

.. _profile-dataframe-export:

Exporting Profiles to DataFrame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One-dimensional profile data can be exported to a :class:`~pandas.DataFrame` object
using the :meth:`yt.data_objects.profiles.Profile1D.to_dataframe` method. Bins which
do not have data will have their fields filled with ``NaN``, except for the bin field
itself. If you only want to export the bins which are used, set ``only_used=True``.

.. code-block:: python

    # Adds all of the data to the DataFrame, but non-used bins are filled with NaNs
    df = profile.to_dataframe()
    # Only adds the used bins to the DataFrame
    df_used = profile.to_dataframe(only_used=True)
    # Only adds the density and temperature fields
    df2 = profile.to_dataframe(fields=[("gas", "density"), ("gas", "temperature")])

The :class:`~pandas.DataFrame` can then analyzed and/or written to disk using pandas
methods. Note that unit information is lost in this export.

.. _profile-astropy-export:

Exporting Profiles to QTable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One-dimensional profile data also can be exported to an AstroPy :class:`~astropy.table.QTable`
object. This table can then be written to disk in a number of formats, such as ASCII text
or FITS files, and manipulated in a number of ways. Bins which do not have data
will have their mask values set to ``False``. If you only want to export the bins
which are used, set ``only_used=True``. Units are preserved in the table by converting
each :class:`~yt.units.yt_array.YTArray` to an :class:`~astropy.units.Quantity`.

To export the 1D profile to a Table object, simply call
:meth:`yt.data_objects.profiles.Profile1D.to_astropy_table`:

.. code-block:: python

    # Adds all of the data to the Table, but non-used bins are masked
    t = profile.to_astropy_table()
    # Only adds the used bins to the Table
    t_used = profile.to_astropy_table(only_used=True)
    # Only adds the density and temperature fields
    t2 = profile.to_astropy_table(fields=[("gas", "density"), ("gas", "temperature")])

.. _generating-line-queries:

Line Queries and Planar Integrals
---------------------------------

To calculate the values along a line connecting two points in a simulation, you
can use the object :class:`~yt.data_objects.selection_data_containers.YTRay`,
accessible as the ``ray`` property on a index.  (See :ref:`data-objects`
for more information on this.)  To do so, you can supply two points and access
fields within the returned object.  For instance, this code will generate a ray
between the points (0.3, 0.5, 0.9) and (0.1, 0.8, 0.5) and examine the density
along that ray:

.. code-block:: python

   ray = ds.ray((0.3, 0.5, 0.9), (0.1, 0.8, 0.5))
   print(ray["gas", "density"])

The points are not ordered, so you may need to sort the data (see the
example in the
:class:`~yt.data_objects.selection_data_containers.YTRay` docs).  Also
note, the ray is traversing cells of varying length, as well as
taking a varying distance to cross each cell.  To determine the
distance traveled by the ray within each cell (for instance, for
integration) the field ``dt`` is available; this field will sum to
1.0, as the ray's path will be normalized to 1.0, independent of how
far it travels through the domain.  To determine the value of ``t`` at
which the ray enters each cell, the field ``t`` is available.  For
instance:

.. code-block:: python

   print(ray["dts"].sum())
   print(ray["t"])

These can be used as inputs to, for instance, the Matplotlib function
:func:`~matplotlib.pyplot.plot`, or they can be saved to disk.

The volume rendering functionality in yt can also be used to calculate
off-axis plane integrals, using the
:class:`~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction`
in a manner similar to that described in :ref:`volume_rendering`.

.. _generating-xarray:

Regular Grids to xarray
-----------------------

Objects that subclass from
:class:`~yt.data_objects.construction_data_containers.YTCoveringGrid` are able
to export to `xarray <https://xarray.pydata.org/>`_.  This enables
interoperability with anything that can take xarray data.  The classes that can do this are
:class:`~yt.data_objects.construction_data_containers.YTCoveringGrid`,
:class:`~yt.data_objects.construction_data_containers.YTArbitraryGrid`, and
:class:`~yt.data_objects.construction_data_containers.YTSmoothedCoveringGrid`.  For example, you can:

.. code-block:: python

   grid = ds.r[::256j, ::256j, ::256j]
   obj = grid.to_xarray(fields=[("gas", "density"), ("gas", "temperature")])

The returned object, ``obj``, will now have the correct labelled axes and so forth.

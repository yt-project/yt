.. _loading-data:

Loading Data
============

This section contains information on how to load data into yt, as well as
some important caveats about different data formats.

.. _loading-art-data:

ART Data
--------

ART data has been supported in the past by Christopher Moody and is currently
cared for by Kenza Arraki.  Please contact the ``yt-dev`` mailing list if you
are interested in using yt for ART data, or if you are interested in assisting
with development of yt to work with ART data.

To load an ART dataset you can use the ``yt.load`` command and provide it the
gas mesh file. It will search for and attempt to find the complementary dark
matter and stellar particle header and data files. However, your simulations may
not follow the same naming convention.

.. code-block:: python

   import yt

   ds = yt.load("D9p_500/10MpcBox_HartGal_csf_a0.500.d")


It will search for and attempt to find the complementary dark matter and stellar
particle header and data files. However, your simulations may not follow the
same naming convention.

For example, the single snapshot given in the sample data has a series of files
that look like this:

.. code-block:: none

   10MpcBox_HartGal_csf_a0.500.d  #Gas mesh
   PMcrda0.500.DAT                #Particle header
   PMcrs0a0.500.DAT               #Particle data (positions,velocities)
   stars_a0.500.dat               #Stellar data (metallicities, ages, etc.)

The ART frontend tries to find the associated files matching the
above, but if that fails you can specify ``file_particle_header``,
``file_particle_data``, and ``file_particle_stars``, in addition to
specifying the gas mesh. Note that the ``pta0.500.dat`` or ``pt.dat``
file containing particle time steps is not loaded by yt.

You also have the option of gridding particles and assigning them onto the
meshes.  This process is in beta, and for the time being it's probably best to
leave ``do_grid_particles=False`` as the default.

To speed up the loading of an ART file, you have a few options. You can turn
off the particles entirely by setting ``discover_particles=False``. You can
also only grid octs up to a certain level, ``limit_level=5``, which is useful
when debugging by artificially creating a 'smaller' dataset to work with.

Finally, when stellar ages are computed we 'spread' the ages evenly within a
smoothing window. By default this is turned on and set to 10Myr. To turn this
off you can set ``spread=False``, and you can tweak the age smoothing window
by specifying the window in seconds, ``spread=1.0e7*365*24*3600``.

There is currently preliminary support for dark matter only ART data. To load a
dataset use the ``yt.load`` command and provide it the particle data file. It
will search for the complementary particle header file.

.. code-block:: python

   import yt

   ds = yt.load("PMcrs0a0.500.DAT")

Important: This should not be used for loading just the dark matter
data for a 'regular' hydrodynamical data set as the units and IO are
different!


.. _loading-artio-data:

ARTIO Data
----------

ARTIO data has a well-specified internal parameter system and has few free
parameters.  However, for optimization purposes, the parameter that provides
the most guidance to yt as to how to manage ARTIO data is ``max_range``.  This
governs the maximum number of space-filling curve cells that will be used in a
single "chunk" of data read from disk.  For small datasets, setting this number
very large will enable more data to be loaded into memory at any given time;
for very large datasets, this parameter can be left alone safely.  By default
it is set to 1024; it can in principle be set as high as the total number of
SFC cells.

To load ARTIO data, you can specify a command such as this:

.. code-block:: python

   ds = load("./A11QR1/s11Qzm1h2_a1.0000.art")

.. _loading_athena_data:

Athena Data
-----------

Athena 4.x VTK data is *mostly* supported and cared for by John
ZuHone. Both uniform grid and SMR datasets are supported.

.. note:
   yt also recognizes Fargo3D data written to VTK files as
   Athena data, but support for Fargo3D data is preliminary.

Loading Athena datasets is slightly different depending on whether
your dataset came from a serial or a parallel run. If the data came
from a serial run or you have joined the VTK files together using the
Athena tool ``join_vtk``, you can load the data like this:

.. code-block:: python

   import yt
   ds = yt.load("kh.0010.vtk")

The filename corresponds to the file on SMR level 0, whereas if there
are multiple levels the corresponding files will be picked up
automatically, assuming they are laid out in ``lev*`` subdirectories
under the directory where the base file is located.

For parallel datasets, yt assumes that they are laid out in
directories named ``id*``, one for each processor number, each with
``lev*`` subdirectories for additional refinement levels. To load this
data, call ``load`` with the base file in the ``id0`` directory:

.. code-block:: python

   import yt
   ds = yt.load("id0/kh.0010.vtk")

which will pick up all of the files in the different ``id*`` directories for
the entire dataset.

yt works in cgs ("Gaussian") units by default, but Athena data is not
normally stored in these units. If you would like to convert data to
cgs units, you may supply conversions for length, time, and mass to ``load`` using
the ``units_override`` functionality:

.. code-block:: python

   import yt

   units_override = {"length_unit":(1.0,"Mpc"),
                     "time_unit"(1.0,"Myr"),
                     "mass_unit":(1.0e14,"Msun")}

   ds = yt.load("id0/cluster_merger.0250.vtk", units_override=units_override)

This means that the yt fields, e.g. ``("gas","density")``, ``("gas","x-velocity")``,
``("gas","magnetic_field_x")``, will be in cgs units, but the Athena fields, e.g.,
``("athena","density")``, ``("athena","velocity_x")``, ``("athena","cell_centered_B_x")``, will be
in code units.

Some 3D Athena outputs may have large grids (especially parallel datasets subsequently joined with
the `join_vtk` script), and may benefit from being subdivided into "virtual grids". For this purpose,
one can pass in the `nprocs` parameter:

.. code-block:: python

   import yt

   ds = yt.load("sloshing.0000.vtk", nprocs=8)

which will subdivide each original grid into `nprocs` grids.

.. note::

    Virtual grids are only supported (and really only necessary) for 3D data.

Alternative values for the following simulation parameters may be specified using a ``parameters``
dict, accepting the following keys:

* ``Gamma``: ratio of specific heats, Type: Float
* ``geometry``: Geometry type, currently accepts ``"cartesian"`` or ``"cylindrical"``
* ``periodicity``: Is the domain periodic? Type: Tuple of boolean values corresponding to each dimension

.. code-block:: python

   import yt

   parameters = {"gamma":4./3., "geometry":"cylindrical", "periodicity":(False,False,False)}

   ds = yt.load("relativistic_jet_0000.vtk", parameters=parameters)

.. rubric:: Caveats

* yt primarily works with primitive variables. If the Athena
  dataset contains conservative variables, the yt primitive fields will be generated from the
  conserved variables on disk.
* Special relativistic datasets may be loaded, but are not fully supported. In particular, the relationships between
  quantities such as pressure and thermal energy will be incorrect, as it is currently assumed that their relationship
  is that of an ideal a :math:`\gamma`-law equation of state.
* Domains may be visualized assuming periodicity.
* Particle list data is currently unsupported.

.. note::

   The old behavior of supplying unit conversions using a ``parameters``
   dict supplied to ``load`` for Athena datasets is still supported, but is being deprecated in
   favor of ``units_override``, which provides the same functionality.

.. _loading-orion-data:

BoxLib Data
-----------

yt has been tested with BoxLib data generated by Orion, Nyx, Maestro and
Castro.  Currently it is cared for by a combination of Andrew Myers, Chris
Malone, Matthew Turk, and Mike Zingale.

To load a BoxLib dataset, you can use the ``yt.load`` command on
the plotfile directory name.  In general, you must also have the
``inputs`` file in the base directory, but Maestro and Castro will get
all the necessary parameter information from the ``job_info`` file in
the plotfile directory.  For instance, if you were in a
directory with the following files:

.. code-block:: none

   inputs
   pltgmlcs5600/
   pltgmlcs5600/Header
   pltgmlcs5600/Level_0
   pltgmlcs5600/Level_0/Cell_H
   pltgmlcs5600/Level_1
   pltgmlcs5600/Level_1/Cell_H
   pltgmlcs5600/Level_2
   pltgmlcs5600/Level_2/Cell_H
   pltgmlcs5600/Level_3
   pltgmlcs5600/Level_3/Cell_H
   pltgmlcs5600/Level_4
   pltgmlcs5600/Level_4/Cell_H

You would feed it the filename ``pltgmlcs5600``:

.. code-block:: python

   import yt
   ds = yt.load("pltgmlcs5600")

For Maestro and Castro, you would not need the ``inputs`` file, and you
would have a ``job_info`` file in the plotfile directory.

.. rubric:: Caveats

* yt does not read the Maestro base state (although you can have Maestro
  map it to a full Cartesian state variable before writing the plotfile
  to get around this).  E-mail the dev list if you need this support.
* yt does not know about particles in Maestro.
* For Maestro, yt aliases either "tfromp" or "tfromh to" ``temperature``
  depending on the value of the ``use_tfromp`` runtime parameter.
* For Maestro, some velocity fields like ``velocity_magnitude`` or
  ``mach_number`` will always use the on-disk value, and not have yt
  derive it, due to the complex interplay of the base state velocity.

.. _loading-pluto-data:

Pluto Data
----------

Support for Pluto AMR data is provided through the Chombo frontend, which
is currently maintained by Andrew Myers. Pluto output files that don't use
the Chombo HDF5 format are currently not supported. To load a Pluto dataset,
you can use the ``yt.load`` command on the ``*.hdf5`` files. For example, the
KelvinHelmholtz sample dataset is a directory that contains the following
files:

.. code-block:: none

   data.0004.hdf5
   pluto.ini

To load it, you can navigate into that directory and do:

.. code-block:: python

   import yt
   ds = yt.load("data.0004.hdf5")

The ``pluto.ini`` file must also be present alongside the HDF5 file.
By default, all of the Pluto fields will be in code units.

.. _loading-enzo-data:

Enzo Data
---------

Enzo data is fully supported and cared for by Matthew Turk.  To load an Enzo
dataset, you can use the ``yt.load`` command and provide it the dataset name.
This would be the name of the output file, and it
contains no extension.  For instance, if you have the following files:

.. code-block:: none

   DD0010/
   DD0010/data0010
   DD0010/data0010.index
   DD0010/data0010.cpu0000
   DD0010/data0010.cpu0001
   DD0010/data0010.cpu0002
   DD0010/data0010.cpu0003

You would feed the ``load`` command the filename ``DD0010/data0010`` as
mentioned.

.. code-block:: python

   import yt
   ds = yt.load("DD0010/data0010")

.. rubric:: Caveats

* There are no major caveats for Enzo usage
* Units should be correct, if you utilize standard unit-setting routines.  yt
  will notify you if it cannot determine the units, although this
  notification will be passive.
* 2D and 1D data are supported, but the extraneous dimensions are set to be
  of length 1.0 in "code length" which may produce strange results for volume
  quantities.

.. _loading-exodusii-data:

Exodus II Data
--------------

.. note::
   To load Exodus II data, you need to have the `netcdf4 <http://unidata.github.io/
   netcdf4-python/>`_ python interface installed.

Exodus II is a file format for Finite Element datasets that is used by the MOOSE
framework for file IO. Support for this format (and for unstructured mesh data in
general) is a new feature as of yt 3.3, so while we aim to fully support it, we
also expect there to be some buggy features at present. Currently, yt can visualize
quads, hexes, triangles, and tetrahedral element types at first order. Additionally,
there is experimental support for the high-order visualization of 20-node hex elements.
Development of more high-order visualization capability is a work in progress.

To load an Exodus II dataset, you can use the ``yt.load`` command on the Exodus II
file:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010", step=0)

Because Exodus II datasets can have multiple steps (which can correspond to time steps,
picard iterations, non-linear solve iterations, etc...), you can also specify a step
argument when you load an Exodus II data that defines the index at which to look when
you read data from the file. Omitting this argument is the same as passing in 0, and
setting ``step=-1`` selects the last time output in the file.

You can access the connectivity information directly by doing:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010", step=-1)
   print(ds.index.meshes[0].connectivity_coords)
   print(ds.index.meshes[0].connectivity_indices)
   print(ds.index.meshes[1].connectivity_coords)
   print(ds.index.meshes[1].connectivity_indices)

This particular dataset has two meshes in it, both of which are made of 8-node hexes.
yt uses a field name convention to access these different meshes in plots and data
objects. To see all the fields found in a particular dataset, you can do:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010")
   print(ds.field_list)

This will give you a list of field names like ``('connect1', 'diffused')`` and
``('connect2', 'convected')``. Here, fields labelled with ``'connect1'`` correspond to the
first mesh, and those with ``'connect2'`` to the second, and so on. To grab the value
of the ``'convected'`` variable at all the nodes in the first mesh, for example, you
would do:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010")
   ad = ds.all_data()  # geometric selection, this just grabs everything
   print(ad['connect1', 'convected'])

In this dataset, ``('connect1', 'convected')`` is nodal field, meaning that the field values
are defined at the vertices of the elements. If we examine the shape of the returned array:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010")
   ad = ds.all_data()
   print(ad['connect1', 'convected'].shape)

we see that this mesh has 12480 8-node hexahedral elements, and that we get 8 field values
for each element. To get the vertex positions at which these field values are defined, we
can do, for instance:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010")
   ad = ds.all_data()
   print(ad['connect1', 'vertex_x'])

If we instead look at an element-centered field, like ``('connect1', 'conv_indicator')``,
we get:

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/out.e-s010")
   ad = ds.all_data()
   print(ad['connect1', 'conv_indicator'].shape)

we instead get only one field value per element.

For information about visualizing unstructured mesh data, including Exodus II datasets,
please see :ref:`unstructured-mesh-slices` and :ref:`unstructured_mesh_rendering`.

Displacement Fields
^^^^^^^^^^^^^^^^^^^

Finite element codes often solve for the displacement of each vertex from its
original position as a node variable, rather than updating the actual vertex
positions with time. For analysis and visualization, it is often useful to turn
these displacements on or off, and to be able to scale them arbitrarily to
emphasize certain features of the solution. To allow this, if ``yt`` detects
displacement fields in an Exodus II dataset (using the convention that they will
be named ``disp_x``, ``disp_y``, etc...), it will optionally add these to
the mesh vertex positions for the purposes of visualization. Displacement fields
can be controlled when a dataset is loaded by passing in an optional dictionary
to the ``yt.load`` command. This feature is turned off by default, meaning that
a dataset loaded as

.. code-block:: python

   import yt
   ds = yt.load("MOOSE_sample_data/mps_out.e")

will not include the displacements in the vertex positions. The displacements can
be turned on separately for each mesh in the file by passing in a a tuple of
(scale, offset) pairs for the meshes you want to enable displacements for.
For example, the following code snippet turns displacements on for the second
mesh, but not the first:

.. code-block:: python

    import yt
    ds = yt.load("MOOSE_sample_data/mps_out.e", step=10,
                 displacements={'connect2': (1.0, [0.0, 0.0, 0.0])})

The displacements can also be scaled by an arbitrary factor before they are
added in to the vertex positions. The following code turns on displacements
for both ``connect1`` and ``connect2``, scaling the former by a factor of 5.0
and the later by a factor of 10.0:

.. code-block:: python

    import yt
    ds = yt.load("MOOSE_sample_data/mps_out.e", step=10,
                 displacements={'connect1': (5.0, [0.0, 0.0, 0.0]),
                                'connect2': (10.0, [0.0, 0.0, 0.0])})

Finally, we can also apply an arbitrary offset to the mesh vertices after
the scale factor is applied. For example, the following code scales all
displacements in the second mesh by a factor of 5.0, and then shifts
each vertex in the mesh by 1.0 unit in the z-direction:

.. code-block:: python

    import yt
    ds = yt.load("MOOSE_sample_data/mps_out.e", step=10,
                  displacements={'connect2': (5.0, [0.0, 0.0, 1.0])})

.. _loading-fits-data:

FITS Data
---------

FITS data is *mostly* supported and cared for by John ZuHone. In order to
read FITS data, `AstroPy <http://www.astropy.org>`_ must be installed. FITS
data cubes can be loaded in the same way by yt as other datasets. yt
can read FITS image files that have the following (case-insensitive) suffixes:

* fits
* fts
* fits.gz
* fts.gz

yt can read two kinds of FITS files: FITS image files and FITS binary table files containing
positions, times, and energies of X-ray events.

.. note::

  AstroPy is necessary due to the requirements of both FITS file reading and
  WCS coordinates. Since new releases of `PyFITS <http://www.stsci
  .edu/institute/software_hardware/pyfits>`_ are to be discontinued, individual
  installations of this package and the `PyWCS <http://stsdas.stsci
  .edu/astrolib/pywcs/>`_ package are not supported.

Though a FITS image is composed of a single array in the FITS file,
upon being loaded into yt it is automatically decomposed into grids:

.. code-block:: python

   import yt
   ds = yt.load("m33_hi.fits")
   ds.print_stats()

.. parsed-literal::

   level  # grids         # cells     # cells^3
   ----------------------------------------------
     0	     512	  981940800       994
   ----------------------------------------------
             512	  981940800

yt will generate its own domain decomposition, but the number of grids can be
set manually by passing the ``nprocs`` parameter to the ``load`` call:

.. code-block:: python

   ds = load("m33_hi.fits", nprocs=1024)

Making the Most of yt for FITS Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yt will load data without WCS information and/or some missing header keywords, but the resulting
field information will necessarily be incomplete. For example, field names may not be descriptive,
and units will not be correct. To get the full use out of yt for FITS files, make sure that for
each image the following header keywords have sensible values:

* ``CDELTx``: The pixel width in along axis ``x``
* ``CRVALx``: The coordinate value at the reference position along axis ``x``
* ``CRPIXx``: The reference pixel along axis ``x``
* ``CTYPEx``: The projection type of axis ``x``
* ``CUNITx``: The units of the coordinate along axis ``x``
* ``BTYPE``: The type of the image
* ``BUNIT``: The units of the image

FITS header keywords can easily be updated using AstroPy. For example,
to set the ``BTYPE`` and ``BUNIT`` keywords:

.. code-block:: python

   import astropy.io.fits as pyfits
   f = pyfits.open("xray_flux_image.fits", mode="update")
   f[0].header["BUNIT"] = "cts/s/pixel"
   f[0].header["BTYPE"] = "flux"
   f.flush()
   f.close()

FITS Coordinates
^^^^^^^^^^^^^^^^

For FITS datasets, the unit of ``code_length`` is always the width of one
pixel. yt will attempt to use the WCS information in the FITS header to
construct information about the coordinate system, and provides support for
the following dataset types:

1. Rectilinear 2D/3D images with length units (e.g., Mpc, AU,
   etc.) defined in the ``CUNITx`` keywords
2. 2D images in some celestial coordinate systems (RA/Dec,
   galactic latitude/longitude, defined in the ``CTYPEx``
   keywords), and X-ray binary table event files
3. 3D images with celestial coordinates and a third axis for another
   quantity, such as velocity, frequency, wavelength, etc.
4. 4D images with the first three axes like Case 3, where the slices
   along the 4th axis are interpreted as different fields.

If your data is of the first case, yt will determine the length units based
on the information in the header. If your data is of the second or third
cases, no length units will be assigned, but the world coordinate information
about the axes will be stored in separate fields. If your data is of the fourth
type, the coordinates of the first three axes will be determined according to
cases 1-3.

.. note::

  Linear length-based coordinates (Case 1 above) are only supported if all dimensions
  have the same value for ``CUNITx``. WCS coordinates are only supported for Cases 2-4.

Fields in FITS Datasets
^^^^^^^^^^^^^^^^^^^^^^^

Multiple fields can be included in a FITS dataset in several different ways.
The first way, and the simplest, is if more than one image HDU is
contained within the same file. The field names will be determined by the
value of ``BTYPE`` in the header, and the field units will be determined by
the value of ``BUNIT``. The second way is if a dataset has a fourth axis,
with each slice along this axis corresponding to a different field. In this
case, the field names will be determined by the value of the ``CTYPE4`` keyword
and the index of the slice. So, for example, if ``BTYPE`` = ``"intensity"`` and
``CTYPE4`` = ``"stokes"``, then the fields will be named
``"intensity_stokes_1"``, ``"intensity_stokes_2"``, and so on.

The third way is if auxiliary files are included along with the main file, like so:

.. code-block:: python

   ds = load("flux.fits", auxiliary_files=["temp.fits","metal.fits"])

The image blocks in each of these files will be loaded as a separate field,
provided they have the same dimensions as the image blocks in the main file.

Additionally, fields corresponding to the WCS coordinates will be generated.
based on the corresponding ``CTYPEx`` keywords. When queried, these fields
will be generated from the pixel coordinates in the file using the WCS
transformations provided by AstroPy.

X-ray event data will be loaded as particle fields in yt, but a grid will be constructed from the
WCS information in the FITS header. There is a helper function, ``setup_counts_fields``,
which may be used to make deposited image fields from the event data for different energy bands
(for an example see :ref:`xray_fits`).

.. note::

  Each FITS image from a single dataset, whether from one file or from one of
  multiple files, must have the same dimensions and WCS information as the
  first image in the primary file. If this is not the case,
  yt will raise a warning and will not load this field.

.. _additional_fits_options:

Additional Options
^^^^^^^^^^^^^^^^^^

The following are additional options that may be passed to the ``load`` command when analyzing
FITS data:

``nan_mask``
""""""""""""

FITS image data may include ``NaNs``. If you wish to mask this data out,
you may supply a ``nan_mask`` parameter, which may either be a
single floating-point number (applies to all fields) or a Python dictionary
containing different mask values for different fields:

.. code-block:: python

   # passing a single float
   ds = load("m33_hi.fits", nan_mask=0.0)

   # passing a dict
   ds = load("m33_hi.fits", nan_mask={"intensity":-1.0,"temperature":0.0})

``suppress_astropy_warnings``
"""""""""""""""""""""""""""""

Generally, AstroPy may generate a lot of warnings about individual FITS
files, many of which you may want to ignore. If you want to see these
warnings, set ``suppress_astropy_warnings = False``.

``z_axis_decomp``
"""""""""""""""""

For some applications, decomposing 3D FITS data into grids that span the x-y plane with short
strides along the z-axis may result in a significant improvement in I/O speed. To enable this feature, set ``z_axis_decomp=True``.

``spectral_factor``
"""""""""""""""""""

Often, the aspect ratio of 3D spectral cubes can be far from unity. Because yt
sets the pixel scale as the ``code_length``, certain visualizations (such as
volume renderings) may look extended or distended in ways that are
undesirable. To adjust the width in ``code_length`` of the spectral axis, set
``spectral_factor`` equal to a constant which gives the desired scaling, or set
it to ``"auto"`` to make the width the same as the largest axis in the sky
plane.

Miscellaneous Tools for Use with FITS Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A number of tools have been prepared for use with FITS data that enhance yt's visualization and
analysis capabilities for this particular type of data. These are included in the ``yt.frontends.fits.misc`` module, and can be imported like so:

.. code-block:: python

  from yt.frontends.fits.misc import setup_counts_fields, PlotWindowWCS, ds9_region

``setup_counts_fields``
"""""""""""""""""""""""

This function can be used to create image fields from X-ray counts data in different energy bands:

.. code-block:: python

  ebounds = [(0.1,2.0),(2.0,5.0)] # Energies are in keV
  setup_counts_fields(ds, ebounds)

which would make two fields, ``"counts_0.1-2.0"`` and ``"counts_2.0-5.0"``,
and add them to the field registry for the dataset ``ds``.


``ds9_region``
""""""""""""""

This function takes a `ds9 <http://ds9.si.edu/site/Home.html>`_ region and creates a "cut region"
data container from it, that can be used to select the cells in the FITS dataset that fall within
the region. To use this functionality, the `pyregion <https://github.com/astropy/pyregion/>`_
package must be installed.

.. code-block:: python

  ds = yt.load("m33_hi.fits")
  circle_region = ds9_region(ds, "circle.reg")
  print(circle_region.quantities.extrema("flux"))


``PlotWindowWCS``
"""""""""""""""""

This class takes a on-axis ``SlicePlot`` or ``ProjectionPlot`` of FITS data and adds celestial
coordinates to the plot axes. To use it, the `WCSAxes <http://wcsaxes.readthedocs.org>`_
package must be installed.

.. code-block:: python

  wcs_slc = PlotWindowWCS(slc)
  wcs_slc.show() # for the IPython notebook
  wcs_slc.save()

``WCSAxes`` is still in an experimental state, but as its functionality improves it will be
utilized more here.

``create_spectral_slabs``
"""""""""""""""""""""""""

.. note::

  The following functionality requires the `spectral-cube <http://spectral-cube.readthedocs.org>`_
  library to be installed.

If you have a spectral intensity dataset of some sort, and would like to extract emission in
particular slabs along the spectral axis of a certain width, ``create_spectral_slabs`` can be
used to generate a dataset with these slabs as different fields. In this example, we use it
to extract individual lines from an intensity cube:

.. code-block:: python

  slab_centers = {'13CN': (218.03117, 'GHz'),
                  'CH3CH2CHO': (218.284256, 'GHz'),
                  'CH3NH2': (218.40956, 'GHz')}
  slab_width = (0.05, "GHz")
  ds = create_spectral_slabs("intensity_cube.fits",
                                    slab_centers, slab_width,
                                    nan_mask=0.0)

All keyword arguments to `create_spectral_slabs` are passed on to `load` when creating the dataset
(see :ref:`additional_fits_options` above). In the returned dataset, the different slabs will be
different fields, with the field names taken from the keys in ``slab_centers``. The WCS coordinates
on the spectral axis are reset so that the center of the domain along this axis is zero, and the
left and right edges of the domain along this axis are :math:`\pm` ``0.5*slab_width``.

Examples of Using FITS Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following IPython notebooks show examples of working with FITS data in yt,
which we recommend you look at in the following order:

* :ref:`radio_cubes`
* :ref:`xray_fits`

.. _loading-flash-data:

FLASH Data
----------

FLASH HDF5 data is *mostly* supported and cared for by John ZuHone.  To load a
FLASH dataset, you can use the ``yt.load`` command and provide it the file name of 
a plot file, checkpoint file, or particle file. Particle files require special handling
depending on the situation, the main issue being that they typically lack grid information. 
The first case is when you have a plotfile and a particle file that you would like to 
load together. In the simplest case, this occurs automatically. For instance, if you
were in a directory with the following files:

.. code-block:: none

   radio_halo_1kpc_hdf5_plt_cnt_0100 # plotfile
   radio_halo_1kpc_hdf5_part_0100 # particle file

where the plotfile and the particle file were created at the same time (therefore having 
particle data consistent with the grid structure of the former). Notice also that the 
prefix ``"radio_halo_1kpc_"`` and the file number ``100`` are the same. In this special case,
the particle file will be loaded automatically when ``yt.load`` is called on the plotfile.
This also works when loading a number of files in a time series.

If the two files do not have the same prefix and number, but they nevertheless have the same
grid structure and are at the same simulation time, the particle data may be loaded with the
``particle_filename`` optional argument to ``yt.load``:

.. code-block:: python

    import yt
    ds = yt.load("radio_halo_1kpc_hdf5_plt_cnt_0100", particle_filename="radio_halo_1kpc_hdf5_part_0100")

However, if you don't have a corresponding plotfile for a particle file, but would still 
like to load the particle data, you can still call ``yt.load`` on the file. However, the 
grid information will not be available, and the particle data will be loaded in a fashion
similar to SPH data. 

.. rubric:: Caveats

* Please be careful that the units are correctly utilized; yt assumes cgs by default, but conversion to
  other :ref:`unit systems <unit_systems>` is also possible.

.. _loading-gadget-data:

Gadget Data
-----------

yt has support for reading Gadget data in both raw binary and HDF5 formats.  It
is able to access the particles as it would any other particle dataset, and it
can apply smoothing kernels to the data to produce both quantitative analysis
and visualization. See :ref:`loading-sph-data` for more details and
:ref:`gadget-notebook` for a detailed example of loading, analyzing, and
visualizing a Gadget dataset.  An example which makes use of a Gadget snapshot
from the OWLS project can be found at :ref:`owls-notebook`.

Gadget data in HDF5 format can be loaded with the ``load`` command:

.. code-block:: python

   import yt
   ds = yt.load("snapshot_061.hdf5")

Gadget data in raw binary format can also be loaded with the ``load`` command.
This is only supported for snapshots created with the ``SnapFormat`` parameter
set to 1 (the standard for Gadget-2).

.. code-block:: python

   import yt
   ds = yt.load("snapshot_061")

.. _particle-bbox:

Units and Bounding Boxes
^^^^^^^^^^^^^^^^^^^^^^^^

There are two additional pieces of information that may be needed.  If your
simulation is cosmological, yt can often guess the bounding box and the units
of the simulation.  However, for isolated simulations and for cosmological
simulations with non-standard units, these must be supplied.  For example, if
a length unit of 1.0 corresponds to a kiloparsec, you can supply this in the
constructor.  yt can accept units such as ``Mpc``, ``kpc``, ``cm``, ``Mpccm/h``
and so on.  In particular, note that ``Mpc/h`` and ``Mpccm/h`` (``cm`` for
comoving here) are usable unit definitions.

yt will attempt to use units for ``mass``, ``length`` and ``time`` as supplied
in the argument ``unit_base``.  The ``bounding_box`` argument is a list of
two-item tuples or lists that describe the left and right extents of the
particles.

.. code-block:: python

   ds = GadgetDataset("snap_004",
           unit_base = {'length': ('kpc', 1.0)},
           bounding_box = [[-600.0, 600.0], [-600.0, 600.0], [-600.0, 600.0]])

.. _particle-indexing-criteria:

Indexing Criteria
^^^^^^^^^^^^^^^^^

yt generates a global mesh index via octree that governs the resolution of
volume elements.  This is governed by two parameters, ``n_ref`` and
``over_refine_factor``.  They are weak proxies for each other.  The first,
``n_ref``, governs how many particles in an oct results in that oct being
refined into eight child octs.  Lower values mean higher resolution; the
default is 64.  The second parameter, ``over_refine_factor``, governs how many
cells are in a given oct; the default value of 1 corresponds to 8 cells.
The number of cells in an oct is defined by the expression
``2**(3*over_refine_factor)``.

It's recommended that if you want higher-resolution, try reducing the value of
``n_ref`` to 32 or 16.

Also yt can be set to generate the global mesh index according to a specific
type of particles instead of all the particles through the parameter
``index_ptype``. For example, to build the octree only according to the
``"PartType0"`` particles, you can do:

.. code-block:: python

   ds = yt.load("snapshot_061.hdf5", index_ptype="PartType0")

By default, ``index_ptype`` is set to ``"all"``, which means all the particles.
Currently this feature only works for the Gadget HDF5 and OWLS datasets. To
bring the feature to other frontends, it's recommended to refer to this
`PR <https://bitbucket.org/yt_analysis/yt/pull-requests/1985/add-particle-type-aware-octree/diff>`_
for implementation details.

.. _gadget-field-spec:

Field Specifications
^^^^^^^^^^^^^^^^^^^^

Binary Gadget outputs often have additional fields or particle types that are
non-standard from the default Gadget distribution format.  These can be
specified in the call to ``GadgetDataset`` by either supplying one of the
sets of field specifications as a string or by supplying a field specification
itself.  As an example, yt has built-in definitions for ``default`` (the
default) and ``agora_unlv``.  Field specifications must be tuples, and must be
of this format:

.. code-block:: python

   default = ( "Coordinates",
               "Velocities",
               "ParticleIDs",
               "Mass",
               ("InternalEnergy", "Gas"),
               ("Density", "Gas"),
               ("SmoothingLength", "Gas"),
   )

This is the default specification used by the Gadget frontend.  It means that
the fields are, in order, Coordinates, Velocities, ParticleIDs, Mass, and the
fields InternalEnergy, Density and SmoothingLength *only* for Gas particles.
So for example, if you have defined a Metallicity field for the particle type
Halo, which comes right after ParticleIDs in the file, you could define it like
this:

.. code-block:: python

   my_field_def = ( "Coordinates",
               "Velocities",
               "ParticleIDs",
               ("Metallicity", "Halo"),
               "Mass",
               ("InternalEnergy", "Gas"),
               ("Density", "Gas"),
               ("SmoothingLength", "Gas"),
   )

To save time, you can utilize the plugins file for yt and use it to add items
to the dictionary where these definitions are stored.  You could do this like
so:

.. code-block:: python

   from yt.frontends.gadget.definitions import gadget_field_specs
   gadget_field_specs["my_field_def"] = my_field_def

Please also feel free to issue a pull request with any new field
specifications, as we're happy to include them in the main distribution!

.. _gadget-ptype-spec:

Particle Type Definitions
^^^^^^^^^^^^^^^^^^^^^^^^^

In some cases, research groups add new particle types or re-order them.  You
can supply alternate particle types by using the keyword ``ptype_spec`` to the
``GadgetDataset`` call.  The default for Gadget binary data is:

.. code-block:: python

   ( "Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry" )

You can specify alternate names, but note that this may cause problems with the
field specification if none of the names match old names.

.. _gadget-header-spec:

Header Specification
^^^^^^^^^^^^^^^^^^^^

If you have modified the header in your Gadget binary file, you can specify an
alternate header specification with the keyword ``header_spec``.  This can
either be a list of strings corresponding to individual header types known to
yt, or it can be a combination of strings and header specifications.  The
default header specification (found in ``yt/frontends/sph/definitions.py``) is:

.. code-block:: python

   default      = (('Npart', 6, 'i'),
                   ('Massarr', 6, 'd'),
                   ('Time', 1, 'd'),
                   ('Redshift', 1, 'd'),
                   ('FlagSfr', 1, 'i'),
                   ('FlagFeedback', 1, 'i'),
                   ('Nall', 6, 'i'),
                   ('FlagCooling', 1, 'i'),
                   ('NumFiles', 1, 'i'),
                   ('BoxSize', 1, 'd'),
                   ('Omega0', 1, 'd'),
                   ('OmegaLambda', 1, 'd'),
                   ('HubbleParam', 1, 'd'),
                   ('FlagAge', 1, 'i'),
                   ('FlagMEtals', 1, 'i'),
                   ('NallHW', 6, 'i'),
                   ('unused', 16, 'i'))

These items will all be accessible inside the object ``ds.parameters``, which
is a dictionary.  You can add combinations of new items, specified in the same
way, or alternately other types of headers.  The other string keys defined are
``pad32``, ``pad64``, ``pad128``, and ``pad256`` each of which corresponds to
an empty padding in bytes.  For example, if you have an additional 256 bytes of
padding at the end, you can specify this with:

.. code-block:: python

   header_spec = ["default", "pad256"]

This can then be supplied to the constructor.  Note that you can also do this
manually, for instance with:


.. code-block:: python

   header_spec = ["default", (('some_value', 8, 'd'),
                              ('another_value', 1, 'i'))]

The letters correspond to data types from the Python struct module.  Please
feel free to submit alternate header types to the main yt repository.

.. _specifying-gadget-units:

Specifying Units
^^^^^^^^^^^^^^^^

If you are running a cosmology simulation, yt will be able to guess the units
with some reliability.  However, if you are not and you do not specify a
dataset, yt will not be able to and will use the defaults of length
being 1.0 Mpc/h (comoving), velocity being in cm/s, and mass being in 10^10
Msun/h.  You can specify alternate units by supplying the ``unit_base`` keyword
argument of this form:

.. code-block:: python

   unit_base = {'length': (1.0, 'cm'), 'mass': (1.0, 'g'), 'time': (1.0, 's')}

yt will utilize length, mass and time to set up all other units.

.. _loading-gamer-data:

GAMER Data
----------

GAMER HDF5 data is supported and cared for by Hsi-Yu Schive. You can load the data like this:

.. code-block:: python

   import yt
   ds = yt.load("InteractingJets/jet_000002")

Currently GAMER does not assume any unit for non-cosmological simulations. To specify the units for yt,
you need to supply conversions for length, time, and mass to ``load`` using the ``units_override`` functionality:

.. code-block:: python

   import yt
   code_units = { "length_unit":(1.0,"kpc"),
                  "time_unit"  :(3.08567758096e+13,"s"),
                  "mass_unit"  :(1.4690033e+36,"g") }
   ds = yt.load("InteractingJets/jet_000002", units_override=code_units)

This means that the yt fields, e.g., ``("gas","density")``, will be in cgs units, but the GAMER fields,
e.g., ``("gamer","Dens")``, will be in code units.

.. rubric:: Caveats

* GAMER data in raw binary format (i.e., OPT__OUTPUT_TOTAL = C-binary) is not supported.

.. _loading-amr-data:

Generic AMR Data
----------------

See :ref:`loading-numpy-array` and
:func:`~yt.frontends.stream.data_structures.load_amr_grids` for more detail.

It is possible to create native yt dataset from Python's dictionary
that describes set of rectangular patches of data of possibly varying
resolution.

.. code-block:: python

   import yt

   grid_data = [
       dict(left_edge = [0.0, 0.0, 0.0],
            right_edge = [1.0, 1.0, 1.],
            level = 0,
            dimensions = [32, 32, 32],
            number_of_particles = 0)
       dict(left_edge = [0.25, 0.25, 0.25],
            right_edge = [0.75, 0.75, 0.75],
            level = 1,
            dimensions = [32, 32, 32],
            number_of_particles = 0)
   ]

   for g in grid_data:
       g["density"] = np.random.random(g["dimensions"]) * 2**g["level"]

   ds = yt.load_amr_grids(grid_data, [32, 32, 32], 1.0)

Particle fields are supported by adding 1-dimensional arrays and
setting the ``number_of_particles`` key to each ``grid``'s dict:

.. code-block:: python

   for g in grid_data:
       g["number_of_particles"] = 100000
       g["particle_position_x"] = np.random.random((g["number_of_particles"]))

.. rubric:: Caveats

* Some functions may behave oddly, and parallelism will be disappointing or
  non-existent in most cases.
* No consistency checks are performed on the index
* Data must already reside in memory.
* Consistency between particle positions and grids is not checked;
  ``load_amr_grids`` assumes that particle positions associated with one grid are
  not bounded within another grid at a higher level, so this must be
  ensured by the user prior to loading the grid data.

Generic Array Data
------------------

See :ref:`loading-numpy-array` and
:func:`~yt.frontends.stream.data_structures.load_uniform_grid` for more detail.

Even if your data is not strictly related to fields commonly used in
astrophysical codes or your code is not supported yet, you can still feed it to
yt to use its advanced visualization and analysis facilities. The only
requirement is that your data can be represented as one or more uniform, three
dimensional numpy arrays. Assuming that you have your data in ``arr``,
the following code:

.. code-block:: python

   import yt

   data = dict(Density = arr)
   bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [1.5, 1.5]])
   ds = yt.load_uniform_grid(data, arr.shape, 3.08e24, bbox=bbox, nprocs=12)

will create yt-native dataset ``ds`` that will treat your array as
density field in cubic domain of 3 Mpc edge size (3 * 3.08e24 cm) and
simultaneously divide the domain into 12 chunks, so that you can take advantage
of the underlying parallelism.

Particle fields are detected as one-dimensional fields. The number of
particles is set by the ``number_of_particles`` key in
``data``. Particle fields are then added as one-dimensional arrays in
a similar manner as the three-dimensional grid fields:

.. code-block:: python

   import yt

   data = dict(Density = dens,
               number_of_particles = 1000000,
               particle_position_x = posx_arr,
	       particle_position_y = posy_arr,
	       particle_position_z = posz_arr)
   bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [1.5, 1.5]])
   ds = yt.load_uniform_grid(data, arr.shape, 3.08e24, bbox=bbox, nprocs=12)

where in this example the particle position fields have been assigned.
``number_of_particles`` must be the same size as the particle arrays. If no
particle arrays are supplied then ``number_of_particles`` is assumed to be
zero.

.. rubric:: Caveats

* Particles may be difficult to integrate.
* Data must already reside in memory.

Semi-Structured Grid Data
-------------------------

See :ref:`loading-numpy-array`,
:func:`~yt.frontends.stream.data_structures.hexahedral_connectivity`,
:func:`~yt.frontends.stream.data_structures.load_hexahedral_mesh` for
more detail.

In addition to uniform grids as described above, you can load in data
with non-uniform spacing between datapoints. To load this type of
data, you must first specify a hexahedral mesh, a mesh of six-sided
cells, on which it will live. You define this by specifying the x,y,
and z locations of the corners of the hexahedral cells. The following
code:

.. code-block:: python

   import yt
   import numpy

   xgrid = numpy.array([-1, -0.65, 0, 0.65, 1])
   ygrid = numpy.array([-1, 0, 1])
   zgrid = numpy.array([-1, -0.447, 0.447, 1])

   coordinates,connectivity = yt.hexahedral_connectivity(xgrid,ygrid,zgrid)

will define the (x,y,z) coordinates of the hexahedral cells and
information about that cell's neighbors such that the cell corners
will be a grid of points constructed as the Cartesion product of
xgrid, ygrid, and zgrid.

Then, to load your data, which should be defined on the interiors of
the hexahedral cells, and thus should have the shape,
``(len(xgrid)-1, len(ygrid)-1, len(zgrid)-1)``, you can use the following code:

.. code-block:: python

   bbox = numpy.array([[numpy.min(xgrid),numpy.max(xgrid)],
	               [numpy.min(ygrid),numpy.max(ygrid)],
	               [numpy.min(zgrid),numpy.max(zgrid)]])
   data = {"density" : arr}
   ds = yt.load_hexahedral_mesh(data,conn,coords,1.0,bbox=bbox)

to load your data into the dataset ``ds`` as described above, where we
have assumed your data is stored in the three-dimensional array
``arr``.

.. rubric:: Caveats

* Integration is not implemented.
* Some functions may behave oddly or not work at all.
* Data must already reside in memory.

Unstructured Grid Data
----------------------

See :ref:`loading-numpy-array`,
:func:`~yt.frontends.stream.data_structures.load_unstructured_mesh` for
more detail.

In addition to the above grid types, you can also load data stored on
unstructured meshes. This type of mesh is used, for example, in many
finite element calculations. Currently, hexahedral and tetrahedral
mesh elements are supported.

To load an unstructured mesh, you need to specify the following. First,
you need to have a coordinates array, which should be an (L, 3) array
that stores the (x, y, z) positions of all of the vertices in the mesh.
Second, you need to specify a connectivity array, which describes how
those vertices are connected into mesh elements. The connectivity array
should be (N, M), where N is the number of elements and M is the
connectivity length, i.e. the number of vertices per element. Finally,
you must also specify a data dictionary, where the keys should be
the names of the fields and the values should be numpy arrays that
contain the field data. These arrays can either supply the cell-averaged
data for each element, in which case they would be (N, 1), or they
can have node-centered data, in which case they would also be (N, M).

Here is an example of how to load an in-memory, unstructured mesh dataset:

.. code-block:: python

   import yt
   import numpy
   from yt.utilities.exodusII_reader import get_data

   coords, connectivity, data = get_data("MOOSE_sample_data/out.e-s010")

This uses a publically available `MOOSE <http://mooseframework.org/>`
dataset along with the get_data function to parse the coords, connectivity,
and data. Then, these can be loaded as an in-memory dataset as follows:

.. code-block:: python

    mesh_id = 0
    ds = yt.load_unstructured_mesh(data[mesh_id], connectivity[mesh_id], coords[mesh_id])

Note that load_unstructured_mesh can take either a single or a list of meshes.
Here, we have selected only the first mesh to load.

.. rubric:: Caveats

* Integration is not implemented.
* Some functions may behave oddly or not work at all.
* Data must already reside in memory.

Generic Particle Data
---------------------

See :ref:`generic-particle-data` and
:func:`~yt.frontends.stream.data_structures.load_particles` for more detail.

You can also load generic particle data using the same ``stream`` functionality
discussed above to load in-memory grid data.  For example, if your particle
positions and masses are stored in ``positions`` and ``massess``, a
vertically-stacked array of particle x,y, and z positions, and a 1D array of
particle masses respectively, you would load them like this:

.. code-block:: python

    import yt

    data = dict(particle_position=positions, particle_mass=masses)
    ds = yt.load_particles(data)

You can also load data using 1D x, y, and z position arrays:

.. code-block:: python

    import yt

    data = dict(particle_position_x=posx,
                particle_position_y=posy,
                particle_position_z=posz,
                particle_mass=masses)
    ds = yt.load_particles(data)

The ``load_particles`` function also accepts the following keyword parameters:

``length_unit``
      The units used for particle positions.

``mass_unit``
       The units of the particle masses.

``time_unit``
       The units used to represent times. This is optional and is only used if
       your data contains a ``creation_time`` field or a ``particle_velocity`` field.

``velocity_unit``
       The units used to represent velocities.  This is optional and is only used
       if you supply a velocity field.  If this is not supplied, it is inferred from
       the length and time units.

``bbox``
       The bounding box for the particle positions.

.. _loading-gizmo-data:

Gizmo Data
----------

Gizmo datasets, including FIRE outputs, can be loaded into yt in the usual 
manner.  Like other SPH data formats, yt loads Gizmo data as particle fields 
and then uses smoothing kernels to deposit those fields to an underlying 
grid structure as spatial fields as described in :ref:`loading-gadget-data`.  
To load Gizmo datasets using the standard HDF5 output format::

   import yt
   ds = yt.load("snapshot_600.hdf5")

Because the Gizmo output format is similar to the Gadget format, yt
may load Gizmo datasets as Gadget depending on the circumstances, but this
should not pose a problem in most situations.  FIRE outputs will be loaded 
accordingly due to the number of metallicity fields found (11 or 17).  

For Gizmo outputs written as raw binary outputs, you may have to specify
a bounding box, field specification, and units as are done for standard 
Gadget outputs.  See :ref:`loading-gadget-data` for more information.

.. _halo-catalog-data:

Halo Catalog Data
-----------------

yt has support for reading halo catalogs produced by Rockstar and the inline
FOF/SUBFIND halo finders of Gadget and OWLS.  The halo catalogs are treated as
particle datasets where each particle represents a single halo.  For example,
this means that the `particle_mass` field refers to the mass of the halos.  For
Gadget FOF/SUBFIND catalogs, the member particles for a given halo can be
accessed by creating `halo` data containers.  See :ref:`halo_containers` for
more information.

If you have access to both the halo catalog and the simulation snapshot from
the same redshift, additional analysis can be performed for each halo using
:ref:`halo_catalog`.

.. _rockstar:

Rockstar
^^^^^^^^

Rockstar halo catalogs are loaded by providing the path to one of the .bin files.
In the case where multiple files were produced, one need only provide the path
to a single one of them.  The field type for all fields is "halos".  Some fields
of note available from Rockstar are:

+----------------+---------------------------+
| Rockstar field | yt field name             |
+================+===========================+
| halo id        | particle_identifier       |
+----------------+---------------------------+
| virial mass    | particle_mass             |
+----------------+---------------------------+
| virial radius  | virial_radius             |
+----------------+---------------------------+
| halo position  | particle_position_(x,y,z) |
+----------------+---------------------------+
| halo velocity  | particle_velocity_(x,y,z) |
+----------------+---------------------------+

Numerous other Rockstar fields exist.  To see them, check the field list by
typing `ds.field_list` for a dataset loaded as `ds`.  Like all other datasets,
fields must be accessed through :ref:`Data-objects`.

.. code-block:: python

   import yt
   ds = yt.load("rockstar_halos/halos_0.0.bin")
   ad = ds.all_data()
   # halo masses
   print(ad["halos", "particle_mass"])
   # halo radii
   print(ad["halos", "virial_radius"])

.. _gadget_fof:

Gadget FOF/SUBFIND
^^^^^^^^^^^^^^^^^^

Gadget FOF/SUBFIND halo catalogs work in the same way as those created by
:ref:`rockstar`, except there are two field types: `FOF` for friend-of-friends
groups and `Subhalo` for halos found with the SUBFIND substructure finder.
Also like Rockstar, there are a number of fields specific to these halo
catalogs.

+-------------------+---------------------------+
| FOF/SUBFIND field | yt field name             |
+===================+===========================+
| halo id           | particle_identifier       |
+-------------------+---------------------------+
| halo mass         | particle_mass             |
+-------------------+---------------------------+
| halo position     | particle_position_(x,y,z) |
+-------------------+---------------------------+
| halo velocity     | particle_velocity_(x,y,z) |
+-------------------+---------------------------+
| num. of particles | particle_number           |
+-------------------+---------------------------+
| num. of subhalos  | subhalo_number (FOF only) |
+-------------------+---------------------------+

Many other fields exist, especially for SUBFIND subhalos.  Check the field
list by typing `ds.field_list` for a dataset loaded as `ds`.  Like all
other datasets, fields must be accessed through :ref:`Data-objects`.

.. code-block:: python

   import yt
   ds = yt.load("gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5")
   ad = ds.all_data()
   # The halo mass
   print(ad["Group", "particle_mass"])
   print(ad["Subhalo", "particle_mass"])
   # Halo ID
   print(ad["Group", "particle_identifier"])
   print(ad["Subhalo", "particle_identifier"])
   # positions
   print(ad["Group", "particle_position_x"])
   # velocities
   print(ad["Group", "particle_velocity_x"])

Multidimensional fields can be accessed through the field name followed by an
underscore and the index.

.. code-block:: python

   # x component of the spin
   print(ad["Subhalo", "SubhaloSpin_0"])

.. _halo_containers:

Halo Data Containers
^^^^^^^^^^^^^^^^^^^^

Halo member particles are accessed by creating halo data containers with the
type of halo ("Group" or "Subhalo") and the halo id.  Scalar values for halos
can be accessed in the same way.  Halos also have mass, position, and velocity
attributes.

.. code-block:: python

   halo = ds.halo("Group", 0)
   # member particles for this halo
   print halo["member_ids"]
   # halo virial radius
   print halo["Group_R_Crit200"]
   # halo mass
   print halo.mass

Subhalos containers can be created using either their absolute ids or their
subhalo ids.

.. code-block:: python

   # first subhalo of the first halo
   subhalo = ds.halo("Subhalo", (0, 0))
   # this subhalo's absolute id
   print subhalo.group_identifier
   # member particles
   print subhalo["member_ids"]

OWLS FOF/SUBFIND
^^^^^^^^^^^^^^^^

OWLS halo catalogs have a very similar structure to regular Gadget halo catalogs.
The two field types are `FOF` and `SUBFIND`.  See :ref:`gadget_fof` for more
information.  At this time, halo member particles cannot be loaded.

.. code-block:: python

   import yt
   ds = yt.load("owls_fof_halos/groups_008/group_008.0.hdf5")
   ad = ds.all_data()
   # The halo mass
   print(ad["FOF", "particle_mass"])

.. _loading-pyne-data:

PyNE Data
---------

`PyNE <http://pyne.io/>`_ is an open source nuclear engineering toolkit
maintained by the PyNE developement team (pyne-dev@googlegroups.com).
PyNE meshes utilize the Mesh-Oriented datABase
`(MOAB) <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB/>`_ and can be
Cartesian or tetrahedral. In addition to field data, pyne meshes store pyne
Material objects which provide a rich set of capabilities for nuclear
engineering tasks. PyNE Cartesian (Hex8) meshes are supported by yt.

To create a pyne mesh:

.. code-block:: python

  from pyne.mesh import Mesh
  num_divisions = 50
  coords = linspace(-1, 1, num_divisions)
  m = Mesh(structured=True, structured_coords=[coords, coords, coords])

Field data can then be added:

.. code-block:: python

  from pyne.mesh import iMeshTag
  m.neutron_flux = IMeshTag()
  # neutron_flux_data is a list or numpy array of size num_divisions^3
  m.neutron_flux[:] = neutron_flux_data

Any field data or material data on the mesh can then be viewed just like any other yt dataset!

.. code-block:: python

  import yt
  pf = yt.frontends.moab.data_structures.PyneMoabHex8Dataset(m)
  s = yt.SlicePlot(pf, 'z', 'neutron_flux')
  s.display()

.. _loading-ramses-data:

RAMSES Data
-----------

In yt-3.0, RAMSES data is fully supported.  If you are interested in taking a
development or stewardship role, please contact the yt-dev mailing list.  To
load a RAMSES dataset, you can use the ``yt.load`` command and provide it
the ``info*.txt`` filename.  For instance, if you were in a
directory with the following files:

.. code-block:: none

   output_00007
   output_00007/amr_00007.out00001
   output_00007/grav_00007.out00001
   output_00007/hydro_00007.out00001
   output_00007/info_00007.txt
   output_00007/part_00007.out00001

You would feed it the filename ``output_00007/info_00007.txt``:

.. code-block:: python

   import yt
   ds = yt.load("output_00007/info_00007.txt")

yt will attempt to guess the fields in the file.  You may also specify a list
of fields by supplying the ``fields`` keyword in your call to ``load``.

.. _loading-sph-data:

SPH Particle Data
-----------------

For all of the SPH frontends, yt uses cython-based SPH smoothing onto an
in-memory octree to create deposited mesh fields from individual SPH particle
fields.

This uses a standard M4 smoothing kernel and the ``smoothing_length``
field to calculate SPH sums, filling in the mesh fields.  This gives you the
ability to both track individual particles (useful for tasks like following
contiguous clouds of gas that would be require a clump finder in grid data) as
well as doing standard grid-based analysis (i.e. slices, projections, and profiles).

The ``smoothing_length`` variable is also useful for determining which particles
can interact with each other, since particles more distant than twice the
smoothing length do not typically see each other in SPH simulations.  By
changing the value of the ``smoothing_length`` and then re-depositing particles
onto the grid, you can also effectively mimic what your data would look like at
lower resolution.

.. _loading-tipsy-data:

Tipsy Data
----------

See :ref:`tipsy-notebook` and :ref:`loading-sph-data` for more details.

yt also supports loading Tipsy data.  Many of its characteristics are similar
to how Gadget data is loaded; specifically, it shares its definition of
indexing and mesh-identification with that described in
:ref:`particle-indexing-criteria`.

.. code-block:: python

   ds = load("./halo1e11_run1.00400")

.. _specifying-cosmology-tipsy:

Specifying Tipsy Cosmological Parameters and Setting Default Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cosmological parameters can be specified to Tipsy to enable computation of
default units.  For example do the following, to load a Tipsy dataset whose
path is stored in the variable ``my_filename`` with specified cosmology
parameters:

.. code-block:: python

   cosmology_parameters = {'current_redshift': 0.0,
                           'omega_lambda': 0.728,
                           'omega_matter': 0.272,
                           'hubble_constant': 0.702}

   ds = yt.load(my_filename,
                cosmology_parameters=cosmology_parameters)

If you wish to set the unit system directly, you can do so by using the
``unit_base`` keyword in the load statement.

 .. code-block:: python

    import yt

    ds = yt.load(filename, unit_base={'length', (1.0, 'Mpc')})

See the documentation for the
:class:`~yt.frontends.tipsy.data_structures.TipsyDataset` class for more
information.

Loading Cosmological Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are not using a parameter file (i.e. non-Gasoline users), then you must
use keyword ``cosmology_parameters`` when loading your data set to indicate to
yt that it is a cosmological data set. If you do not wish to set any
non-default cosmological parameters, you may pass an empty dictionary.

 .. code-block:: python

    import yt
    ds = yt.load(filename, cosmology_parameters={})




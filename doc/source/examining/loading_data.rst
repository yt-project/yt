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
meshes.  This process is in beta, and for the time being, it's probably best to
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

.. _loading-athena-data:

Athena Data
-----------

Athena 4.x VTK data is supported and cared for by John ZuHone. Both uniform grid
and SMR datasets are supported.

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

The default unit system in yt is cgs ("Gaussian") units, but Athena data is not
normally stored in these units, so the code unit system is the default unit
system for Athena data. This means that answers to field queries from data
objects and plots of data will be expressed in code units. Note that the default
conversions from these units will still be in terms of cgs units, e.g. 1
``code_length`` equals 1 cm, and so on. If you would like to provided different
conversions, you may supply conversions for length, time, and mass to ``load``
using the ``units_override`` functionality:

.. code-block:: python

   import yt

   units_override = {"length_unit": (1.0, "Mpc"),
                     "time_unit": (1.0, "Myr"),
                     "mass_unit": (1.0e14, "Msun")}

   ds = yt.load("id0/cluster_merger.0250.vtk", units_override=units_override)

This means that the yt fields, e.g. ``("gas","density")``,
``("gas","velocity_x")``, ``("gas","magnetic_field_x")``, will be in cgs units
(or whatever unit system was specified), but the Athena fields, e.g.,
``("athena","density")``, ``("athena","velocity_x")``,
``("athena","cell_centered_B_x")``, will be in code units.

Some 3D Athena outputs may have large grids (especially parallel datasets
subsequently joined with the ``join_vtk`` script), and may benefit from being
subdivided into "virtual grids". For this purpose, one can pass in the
``nprocs`` parameter:

.. code-block:: python

   import yt

   ds = yt.load("sloshing.0000.vtk", nprocs=8)

which will subdivide each original grid into ``nprocs`` grids. Note that this
parameter is independent of the number of MPI tasks assigned to analyze the data
set in parallel (see :ref:`parallel-computation`), and ideally should be (much)
larger than this.

.. note::

    Virtual grids are only supported (and really only necessary) for 3D data.

Alternative values for the following simulation parameters may be specified
using a ``parameters`` dict, accepting the following keys:

* ``Gamma``: ratio of specific heats, Type: Float
* ``geometry``: Geometry type, currently accepts ``"cartesian"`` or
  ``"cylindrical"``
* ``periodicity``: Is the domain periodic? Type: Tuple of boolean values
  corresponding to each dimension

.. code-block:: python

   import yt

   parameters = {"gamma":4./3., "geometry":"cylindrical",
                 "periodicity":(False,False,False)}

   ds = yt.load("relativistic_jet_0000.vtk", parameters=parameters)

.. rubric:: Caveats

* yt primarily works with primitive variables. If the Athena dataset contains
  conservative variables, the yt primitive fields will be generated from the
  conserved variables on disk.
* Special relativistic datasets may be loaded, but at this time not all of
  their fields are fully supported. In particular, the relationships between
  quantities such as pressure and thermal energy will be incorrect, as it is
  currently assumed that their relationship is that of an ideal a
  :math:`\gamma`-law equation of state. This will be rectified in a future
  release.
* Domains may be visualized assuming periodicity.
* Particle list data is currently unsupported.

.. note::

   The old behavior of supplying unit conversions using a ``parameters``
   dict supplied to ``load`` for Athena datasets is still supported, but is
   being deprecated in favor of ``units_override``, which provides the same
   functionality.

.. _loading-athena-pp-data:

Athena++ Data
-------------

Athena++ HDF5 data is supported and cared for by John ZuHone. Uniform-grid, SMR,
and AMR datasets in cartesian coordinates are fully supported. Support for
curvilinear coordinates and logarithmic cell sizes exists, but is preliminary.
For the latter type of dataset, the data will be loaded in as a semi-structured
mesh dataset. See :ref:`loading-semi-structured-mesh-data` for more details on
how this works in yt.

The default unit system in yt is cgs ("Gaussian") units, but Athena++ data is
not normally stored in these units, so the code unit system is the default unit
system for Athena++ data. This means that answers to field queries from data
objects and plots of data will be expressed in code units. Note that the default
conversions from these units will still be in terms of cgs units, e.g. 1
``code_length`` equals 1 cm, and so on. If you would like to provided different
conversions, you may supply conversions for length, time, and mass to ``load``
using the ``units_override`` functionality:

.. code-block:: python

   import yt

   units_override = {"length_unit":(1.0,"Mpc"),
                     "time_unit"(1.0,"Myr"),
                     "mass_unit":(1.0e14,"Msun")}

   ds = yt.load("AM06/AM06.out1.00400.athdf", units_override=units_override)

This means that the yt fields, e.g. ``("gas","density")``,
``("gas","velocity_x")``, ``("gas","magnetic_field_x")``, will be in cgs units
(or whatever unit system was specified), but the Athena fields, e.g.,
``("athena_pp","density")``, ``("athena_pp","vel1")``, ``("athena_pp","Bcc1")``,
will be in code units.

.. rubric:: Caveats

* yt primarily works with primitive variables. If the Athena++ dataset contains
  conservative variables, the yt primitive fields will be generated from the
  conserved variables on disk.
* Special relativistic datasets may be loaded, but at this time not all of their
  fields are fully supported. In particular, the relationships between
  quantities such as pressure and thermal energy will be incorrect, as it is
  currently assumed that their relationship is that of an ideal
  :math:`\gamma`-law equation of state. This will be rectified in a future
  release.
* Domains may be visualized assuming periodicity.

.. _loading-orion-data:

AMReX / BoxLib Data
-------------------

AMReX and BoxLib share a frontend (currently named `boxlib`), since
the file format nearly identical.  yt has been tested with AMReX/BoxLib
data generated by Orion, Nyx, Maestro, Castro, IAMR, and
WarpX. Currently it is cared for by a combination of Andrew Myers,
Matthew Turk, and Mike Zingale.

To load an AMReX/BoxLib dataset, you can use the ``yt.load`` command on
the plotfile directory name.  In general, you must also have the
``inputs`` file in the base directory, but Maestro, Castro, Nyx, and WarpX will get
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

For Maestro, Castro, Nyx, and WarpX, you would not need the ``inputs`` file, and you
would have a ``job_info`` file in the plotfile directory.

.. rubric:: Caveats

* yt does not read the Maestro base state (although you can have Maestro
  map it to a full Cartesian state variable before writing the plotfile
  to get around this).  E-mail the dev list if you need this support.
* yt supports AMReX/BoxLib particle data stored in the standard format used
  by Nyx and WarpX, and optionally Castro. It currently does not support the ASCII particle
  data used by Maestro and Castro.
* For Maestro, yt aliases either "tfromp" or "tfromh to" ``temperature``
  depending on the value of the ``use_tfromp`` runtime parameter.
* For Maestro, some velocity fields like ``velocity_magnitude`` or
  ``mach_number`` will always use the on-disk value, and not have yt
  derive it, due to the complex interplay of the base state velocity.

Viewing raw fields in WarpX
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most AMReX/BoxLib codes output cell-centered data. If the underlying discretization
is not cell-centered, then fields are typically averaged to cell centers before
they are written to plot files for visualization. WarpX, however, has the option
to output the raw (i.e., not averaged to cell centers) data as well.  If you
run your WarpX simulation with ``warpx.plot_raw_fields = 1`` in your inputs
file, then you should get an additional ``raw_fields`` subdirectory inside your
plot file. When you load this dataset, yt will have additional on-disk fields
defined, with the "raw" field type:

.. code-block:: python

    import yt
    ds = yt.load("Laser/plt00015/")
    print(ds.field_list)

The raw fields in WarpX are nodal in at least one direction. We define a field
to be "nodal" in a given direction if the field data is defined at the "low"
and "high" sides of the cell in that direction, rather than at the cell center.
Instead of returning one field value per cell selected, nodal fields return a
number of values, depending on their centering. This centering is marked by
a `nodal_flag` that describes whether the fields is nodal in each dimension.
``nodal_flag = [0, 0, 0]`` means that the field is cell-centered, while
``nodal_flag = [0, 0, 1]`` means that the field is nodal in the z direction
and cell centered in the others, i.e. it is defined on the z faces of each cell.
``nodal_flag = [1, 1, 0]`` would mean that the field is centered in the z direction,
but nodal in the other two, i.e. it lives on the four cell edges that are normal
to the z direction.

.. code-block:: python

    ds.index
    ad = ds.all_data()
    print(ds.field_info[('raw', 'Ex')].nodal_flag)
    print(ad['raw', 'Ex'].shape)
    print(ds.field_info[('raw', 'Bx')].nodal_flag)
    print(ad['raw', 'Bx'].shape)
    print(ds.field_info[('boxlib', 'Bx')].nodal_flag)
    print(ad['boxlib', 'Bx'].shape)

Here, the field ``('raw', 'Ex')`` is nodal in two directions, so four values per cell
are returned, corresponding to the four edges in each cell on which the variable
is defined. ``('raw', 'Bx')`` is nodal in one direction, so two values are returned
per cell. The standard, averaged-to-cell-centers fields are still available.

Currently, slices and data selection are implemented for nodal fields. Projections,
volume rendering, and many of the analysis modules will not work.

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


Enzo MHDCT data
^^^^^^^^^^^^^^^

The electric and magnetic fields for Enzo MHDCT simulations are defined on cell
faces, unlike other Enzo fields which are defined at cell centers. In yt, we
call face-centered fields like this "nodal".  We define a field to be nodal in
a given direction if the field data is defined at the "low" and "high" sides of
the cell in that direction, rather than at the cell center.  Instead of
returning one field value per cell selected, nodal fields return a number of
values, depending on their centering. This centering is marked by a `nodal_flag`
that describes whether the fields is nodal in each dimension.  ``nodal_flag =
[0, 0, 0]`` means that the field is cell-centered, while ``nodal_flag = [0, 0,
1]`` means that the field is nodal in the z direction and cell centered in the
others, i.e. it is defined on the z faces of each cell.  ``nodal_flag = [1, 1,
0]`` would mean that the field is centered in the z direction, but nodal in the
other two, i.e. it lives on the four cell edges that are normal to the z
direction.

.. code-block:: python

    ds.index
    ad = ds.all_data()
    print(ds.field_info[('enzo', 'Ex')].nodal_flag)
    print(ad['raw', 'Ex'].shape)
    print(ds.field_info[('enzo', 'BxF')].nodal_flag)
    print(ad['raw', 'Bx'].shape)
    print(ds.field_info[('enzo', 'Bx')].nodal_flag)
    print(ad['boxlib', 'Bx'].shape)

Here, the field ``('enzo', 'Ex')`` is nodal in two directions, so four values
per cell are returned, corresponding to the four edges in each cell on which the
variable is defined. ``('enzo', 'BxF')`` is nodal in one direction, so two
values are returned per cell. The standard, non-nodal field ``('enzo', 'Bx')``
is also available.

Currently, slices and data selection are implemented for nodal
fields. Projections, volume rendering, and many of the analysis modules will not
work.

.. _loading-enzop-data:

Enzo-P Data
-----------

Enzo-P outputs have three types of files.

.. code-block:: none

   hello-0200/
   hello-0200/hello-0200.block_list
   hello-0200/hello-0200.file_list
   hello-0200/hello-0200.hello-c0020-p0000.h5

To load Enzo-P data into yt, provide the block list file:

.. code-block:: python

   import yt
   ds = yt.load("hello-0200/hello-0200.block_list")

Mesh and particle fields are fully supported for 1, 2, and 3D datasets.  Enzo-P
supports arbitrary particle types defined by the user.  The available particle
types will be known as soon as the dataset index is created.

.. code-block:: python

   ds = yt.load("ENZOP_DD0140/ENZOP_DD0140.block_list")
   ds.index
   print(ds.particle_types)
   print(ds.particle_type_counts)
   print(ds.r["dark", "particle_position"])

.. rubric:: Caveats

* The Enzo-P output format is still evolving somewhat as the code is being
  actively developed. This frontend will be updated as development continues
  and backward compatibility may occasionally be broken until the file format
  has converged.

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

yt can currently read two kinds of FITS files: FITS image files and FITS
binary table files containing positions, times, and energies of X-ray events.

Though a FITS image is composed of a single array in the FITS file,
upon being loaded into yt it is automatically decomposed into grids:

.. code-block:: python

   import yt
   ds = yt.load("m33_hi.fits")
   ds.print_stats()

.. parsed-literal::

   level  # grids         # cells     # cells^3
   ----------------------------------------------
     0       512          981940800       994
   ----------------------------------------------
             512          981940800

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
about the axes will be stored in separate fields. If your data is of the
fourth type, the coordinates of the first three axes will be determined
according to cases 1-3.

.. note::

  Linear length-based coordinates (Case 1 above) are only supported if all
  dimensions have the same value for ``CUNITx``. WCS coordinates are only
  supported for Cases 2-4.

FITS Data Decomposition
^^^^^^^^^^^^^^^^^^^^^^^

Though a FITS image is composed of a single array in the FITS file,
upon being loaded into yt it is automatically decomposed into grids:

.. code-block:: python

   import yt
   ds = yt.load("m33_hi.fits")
   ds.print_stats()

.. parsed-literal::

   level  # grids         # cells     # cells^3
   ----------------------------------------------
     0       512          981940800       994
   ----------------------------------------------
             512          981940800

For 3D spectral-cube data, the decomposition into grids will be done along the
spectral axis since this will speed up many common operations for this
particular type of dataset.

yt will generate its own domain decomposition, but the number of grids can be
set manually by passing the ``nprocs`` parameter to the ``load`` call:

.. code-block:: python

   ds = load("m33_hi.fits", nprocs=64)


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

X-ray event data will be loaded as particle fields in yt, but a grid will be
constructed from the WCS information in the FITS header. There is a helper
function, ``setup_counts_fields``, which may be used to make deposited image
fields from the event data for different energy bands (for an example see
:ref:`xray_fits`).

.. note::

  Each FITS image from a single dataset, whether from one file or from one of
  multiple files, must have the same dimensions and WCS information as the
  first image in the primary file. If this is not the case,
  yt will raise a warning and will not load this field.

.. _additional_fits_options:

Additional Options
^^^^^^^^^^^^^^^^^^

The following are additional options that may be passed to the ``load`` command
when analyzing FITS data:

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

A number of tools have been prepared for use with FITS data that enhance yt's
visualization and analysis capabilities for this particular type of data. These
are included in the ``yt.frontends.fits.misc`` module, and can be imported like
so:

.. code-block:: python

  from yt.frontends.fits.misc import setup_counts_fields, PlotWindowWCS, ds9_region

``setup_counts_fields``
"""""""""""""""""""""""

This function can be used to create image fields from X-ray counts data in
different energy bands:

.. code-block:: python

  ebounds = [(0.1,2.0),(2.0,5.0)] # Energies are in keV
  setup_counts_fields(ds, ebounds)

which would make two fields, ``"counts_0.1-2.0"`` and ``"counts_2.0-5.0"``,
and add them to the field registry for the dataset ``ds``.

``ds9_region``
""""""""""""""

This function takes a `ds9 <http://ds9.si.edu/site/Home.html>`_ region and
creates a "cut region" data container from it, that can be used to select
the cells in the FITS dataset that fall within the region. To use this
functionality, the `pyregion <https://github.com/astropy/pyregion/>`_
package must be installed.

.. code-block:: python

  ds = yt.load("m33_hi.fits")
  circle_region = ds9_region(ds, "circle.reg")
  print(circle_region.quantities.extrema("flux"))


``PlotWindowWCS``
"""""""""""""""""

This class takes a on-axis ``SlicePlot`` or ``ProjectionPlot`` of FITS
data and adds celestial coordinates to the plot axes. To use it, a
version of AstroPy >= 1.3 must be installed.

.. code-block:: python

  wcs_slc = PlotWindowWCS(slc)
  wcs_slc.show() # for the IPython notebook
  wcs_slc.save()

``WCSAxes`` is still in an experimental state, but as its functionality
improves it will be utilized more here.

``create_spectral_slabs``
"""""""""""""""""""""""""

.. note::

  The following functionality requires the
  `spectral-cube <http://spectral-cube.readthedocs.org>`_ library to be
  installed.

If you have a spectral intensity dataset of some sort, and would like to
extract emission in particular slabs along the spectral axis of a certain
width, ``create_spectral_slabs`` can be used to generate a dataset with
these slabs as different fields. In this example, we use it to extract
individual lines from an intensity cube:

.. code-block:: python

  slab_centers = {'13CN': (218.03117, 'GHz'),
                  'CH3CH2CHO': (218.284256, 'GHz'),
                  'CH3NH2': (218.40956, 'GHz')}
  slab_width = (0.05, "GHz")
  ds = create_spectral_slabs("intensity_cube.fits",
                                    slab_centers, slab_width,
                                    nan_mask=0.0)

All keyword arguments to ``create_spectral_slabs`` are passed on to ``load`` when
creating the dataset (see :ref:`additional_fits_options` above). In the
returned dataset, the different slabs will be different fields, with the field
names taken from the keys in ``slab_centers``. The WCS coordinates on the
spectral axis are reset so that the center of the domain along this axis is
zero, and the left and right edges of the domain along this axis are
:math:`\pm` ``0.5*slab_width``.

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

.. note:: 

   If you are loading a multi-file dataset with Gadget, supply the *zeroth*
   file to the ``load`` command.  For instance,
   ``yt.load("snapshot_061.0.hdf5")`` .

Gadget data in HDF5 format can be loaded with the ``load`` command:

.. code-block:: python

   import yt
   ds = yt.load("snapshot_061.hdf5")

Gadget data in raw binary format can also be loaded with the ``load`` command.
This is supported for snapshots created with the ``SnapFormat`` parameter
set to 1 or 2.

.. code-block:: python

   import yt
   ds = yt.load("snapshot_061")

.. _particle-bbox:

Units and Bounding Boxes
^^^^^^^^^^^^^^^^^^^^^^^^

There are two additional pieces of information that may be needed.  If your
simulation is cosmological, yt can often guess the bounding box and the units of
the simulation.  However, for isolated simulations and for cosmological
simulations with non-standard units, these must be supplied by the user.  For
example, if a length unit of 1.0 corresponds to a kiloparsec, you can supply
this in the constructor.  yt can accept units such as ``Mpc``, ``kpc``, ``cm``,
``Mpccm/h`` and so on.  In particular, note that ``Mpc/h`` and ``Mpccm/h``
(``cm`` for comoving here) are usable unit definitions.

yt will attempt to use units for ``mass``, ``length`` and ``time`` as supplied
in the argument ``unit_base``.  The ``bounding_box`` argument is a list of
two-item tuples or lists that describe the left and right extents of the
particles. In this example we load a dataset with a custom bounding box
and units.

.. code-block:: python


   bbox = [[-600.0, 600.0], [-600.0, 600.0], [-600.0, 600.0]]
   unit_base = {
       'length': (1.0, 'kpc'),
       'velocity: (1.0, 'km/s'),
       'mass': (1.0, 'Msun')
   }

   ds = yt.load("snap_004", unit_base=unit_base, bounding_box=bbox)

In addition, you can use ``UnitLength_in_cm``, ``UnitVelocity_in_cm_per_s``,
and ``UnitMass_in_g`` as keys for the ``unit_base`` dictionary. These names
come from the names used in the Gadget runtime parameter file. This example
will initialize a dataset with the same units as the example above:

.. code-block:: python

  unit_base = {
      'UnitLength_in_cm': 3.09e21,
      'UnitVelocity_in_cm_per_s': 1e5
      'UnitMass_in_g': 1.989e33
   }

  ds = yt.load("snap_004", unit_base=unit_base, bounding_box=bbox)

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
For Gadget binary outputs, ``index_ptype`` should be set using the particle type
names yt uses internally (e.g. ``'Gas'``, ``'Halo'``, ``'Disk'``, etc). For
Gadget HDF5 outputs the particle type names come from the HDF5 output and so
should be referred to using names like ``'PartType0'``.

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

   header_spec = "default+pad256"

Note that a single string like this means a single header block.  To specify
multiple header blocks, use a list of strings instead:

.. code-block:: python

  header_spec = ["default", "pad256"]

This can then be supplied to the constructor.  Note that you can also define
header items manually, for instance with:

.. code-block:: python

   from yt.frontends.gadget.definitions import gadget_header_specs

   gadget_header_specs["custom"] = (('some_value', 8, 'd'),
                                    ('another_value', 1, 'i'))
   header_spec = "default+custom"

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

For simulations without units (i.e., OPT__UNIT = 0), you can supply conversions for
length, time, and mass to ``load`` using the ``units_override`` functionality:

.. code-block:: python

   import yt
   code_units = { "length_unit":(1.0,"kpc"),
                  "time_unit"  :(3.08567758096e+13,"s"),
                  "mass_unit"  :(1.4690033e+36,"g") }
   ds = yt.load("InteractingJets/jet_000002", units_override=code_units)

This means that the yt fields, e.g., ``("gas","density")``, will be in cgs units, but the GAMER fields,
e.g., ``("gamer","Dens")``, will be in code units.

Particle data are supported and are always stored in the same file as the grid data.

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
       dict(left_edge=[0.0, 0.0, 0.0],
            right_edge=[1.0, 1.0, 1.0],
            level=0,
            dimensions=[32, 32, 32])
       dict(left_edge=[0.25, 0.25, 0.25],
            right_edge=[0.75, 0.75, 0.75],
            level=1,
            dimensions=[32, 32, 32])
   ]

   for g in grid_data:
       g["density"] = np.random.random(g["dimensions"]) * 2 ** g["level"]

   ds = yt.load_amr_grids(grid_data, [32, 32, 32], 1.0)

.. note::

   yt only supports a block structure where the grid edges on the ``n``-th
   refinement level are aligned with the cell edges on the ``n-1``-th level.

Particle fields are supported by adding 1-dimensional arrays to each
``grid``'s dict:

.. code-block:: python

   for g in grid_data:
       g["particle_position_x"] = np.random.random(size=100000)

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

Particle fields are added as one-dimensional arrays in a similar manner as the
three-dimensional grid fields:

.. code-block:: python

   import yt

   data = dict(Density = dens,
               particle_position_x = posx_arr,
                   particle_position_y = posy_arr,
                   particle_position_z = posz_arr)
   bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [1.5, 1.5]])
   ds = yt.load_uniform_grid(data, arr.shape, 3.08e24, bbox=bbox, nprocs=12)

where in this example the particle position fields have been assigned. If no
particle fields are supplied, then the number of particles is assumed to be
zero.

.. rubric:: Caveats

* Particles may be difficult to integrate.
* Data must already reside in memory.

.. _loading-semi-structured-mesh-data:

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
will be a grid of points constructed as the Cartesian product of
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
    import numpy as np

    coords = np.array([[0.0, 0.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [0.0, 1.0]], dtype=np.float64)

     connect = np.array([[0, 1, 3],
                         [1, 2, 3]], dtype=np.int64)

     data = {}
     data['connect1', 'test'] = np.array([[0.0, 1.0, 3.0],
                                          [1.0, 2.0, 3.0]], dtype=np.float64)

Here, we have made up a simple, 2D unstructured mesh dataset consisting of two
triangles and one node-centered data field. This data can be loaded as an in-memory
dataset as follows:

.. code-block:: python

    ds = yt.load_unstructured_mesh(connect, coords, data)

The in-memory dataset can then be visualized as usual, e.g.:

.. code-block:: python

    sl = yt.SlicePlot(ds, 'z', 'test')
    sl.annotate_mesh_lines()

Note that load_unstructured_mesh can take either a single mesh or a list of meshes.
To load multiple meshes, you can do:

.. code-block:: python

   import yt
   import numpy as np

   coordsMulti = np.array([[0.0, 0.0],
                           [1.0, 0.0],
                           [1.0, 1.0],
                           [0.0, 1.0]], dtype=np.float64)

   connect1 = np.array([[0, 1, 3], ], dtype=np.int64)
   connect2 = np.array([[1, 2, 3], ], dtype=np.int64)

   data1 = {}
   data2 = {}
   data1['connect1', 'test'] = np.array([[0.0, 1.0, 3.0], ], dtype=np.float64)
   data2['connect2', 'test'] = np.array([[1.0, 2.0, 3.0], ], dtype=np.float64)

   connectList = [connect1, connect2]
   dataList    = [data1, data2]

   ds = yt.load_unstructured_mesh(connectList, coordsMulti, dataList)

   # only plot the first mesh
   sl = yt.SlicePlot(ds, 'z', ('connect1', 'test'))

   # only plot the second
   sl = yt.SlicePlot(ds, 'z', ('connect2', 'test'))

   # plot both
   sl = yt.SlicePlot(ds, 'z', ('all', 'test'))

Note that you must respect the field naming convention that fields on the first
mesh will have the type 'connect1', fields on the second will have 'connect2', etc...

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
positions and masses are stored in ``positions`` and ``masses``, a
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

If ``("PartType0", "MagneticField")`` is present in the output, it would be
loaded and aliased to ``("PartType0", "particle_magnetic_field")``. The
corresponding component field like ``("PartType0", "particle_magnetic_field_x")``
would be added automatically.

Note that ``("PartType4", "StellarFormationTime")`` field has different
meanings depending on whether it is a cosmological simulation. For cosmological
runs this is the scale factor at the redshift when the star particle formed.
For non-cosmological runs it is the time when the star particle formed. (See the
`GIZMO User Guide <http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html>`_)
For this reason, ``("PartType4", "StellarFormationTime")`` is loaded as a
dimensionless field. We defined two related fields
``("PartType4", "creation_time")``, and ``("PartType4", "age")`` with physical
units for your convenience.

For Gizmo outputs written as raw binary outputs, you may have to specify
a bounding box, field specification, and units as are done for standard
Gadget outputs.  See :ref:`loading-gadget-data` for more information.

.. _halo-catalog-data:

Halo Catalog Data
-----------------

yt has support for reading halo catalogs produced by the Amiga Halo Finder (AHF), Rockstar and the inline
FOF/SUBFIND halo finders of Gadget and OWLS.  The halo catalogs are treated as
particle datasets where each particle represents a single halo.  For example,
this means that the `particle_mass` field refers to the mass of the halos.  For
Gadget FOF/SUBFIND catalogs, the member particles for a given halo can be
accessed by creating `halo` data containers.  See :ref:`halo_containers` for
more information.

If you have access to both the halo catalog and the simulation snapshot from
the same redshift, additional analysis can be performed for each halo using
:ref:`halo_catalog`.  The resulting product can be reloaded in a similar manner
to the other halo catalogs shown here.

.. _ahf:

Amiga Halo Finder
^^^^^^^^^^^^^^^^^

Amiga Halo Finder (AHF) halo catalogs are loaded by providing the path to the
.parameter files.  The corresponding .log and .AHF_halos files must exist for
data loading to succeed. The field type for all fields is "halos". Some fields
of note available from AHF are:

+----------------+---------------------------+
| AHF field      | yt field name             |
+================+===========================+
| ID             | particle_identifier       |
+----------------+---------------------------+
| Mvir           | particle_mass             |
+----------------+---------------------------+
| Rvir           | virial_radius             |
+----------------+---------------------------+
| (X,Y,Z)c       | particle_position_(x,y,z) |
+----------------+---------------------------+
| V(X,Y,Z)c      | particle_velocity_(x,y,z) |
+----------------+---------------------------+

Numerous other AHF fields exist.  To see them, check the field list by typing
`ds.field_list` for a dataset loaded as `ds`.  Like all other datasets, fields
must be accessed through :ref:`Data-objects`.

.. code-block:: python

   import yt
   ds = yt.load("ahf_halos/snap_N64L16_135.parameter", hubble_constant=0.7)
   ad = ds.all_data()
   # halo masses
   print(ad["halos", "particle_mass"])
   # halo radii
   print(ad["halos", "virial_radius"])

.. note::

  Currently the dimensionless Hubble parameter that yt needs is not provided in
  AHF outputs. So users need to provide the `hubble_constant` (default to 1.0) while loading datasets, as shown above.

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
   print(halo["member_ids"])
   # halo virial radius
   print(halo["Group_R_Crit200"])
   # halo mass
   print(halo.mass)

Subhalos containers can be created using either their absolute ids or their
subhalo ids.

.. code-block:: python

   # first subhalo of the first halo
   subhalo = ds.halo("Subhalo", (0, 0))
   # this subhalo's absolute id
   print(subhalo.group_identifier)
   # member particles
   print(subhalo["member_ids"])

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

.. _halocatalog:

HaloCatalog
^^^^^^^^^^^

These are catalogs produced by the analysis discussed in :ref:`halo_catalog`.
In the case where multiple files were produced, one need only provide the path
to a single one of them.  The field type for all fields is "halos".  The fields
available here are similar to other catalogs.  Any addition
:ref:`halo_catalog_quantities` will also be accessible as fields.

+-------------------+---------------------------+
| HaloCatalog field | yt field name             |
+===================+===========================+
| halo id           | particle_identifier       |
+-------------------+---------------------------+
| virial mass       | particle_mass             |
+-------------------+---------------------------+
| virial radius     | virial_radius             |
+-------------------+---------------------------+
| halo position     | particle_position_(x,y,z) |
+-------------------+---------------------------+
| halo velocity     | particle_velocity_(x,y,z) |
+-------------------+---------------------------+

.. code-block:: python

   import yt
   ds = yt.load("catalogs/catalog.0.h5")
   ad = ds.all_data()
   # The halo mass
   print(ad["halos", "particle_mass"])

.. _loading-openpmd-data:

openPMD Data
------------

`openPMD <http://www.openpmd.org>`_ is an open source meta-standard and naming
scheme for mesh based data and particle data. It does not actually define a file
format.

HDF5-containers respecting the minimal set of meta information from
versions 1.0.0 and 1.0.1 of the standard are compatible.
Support for the ED-PIC extension is not available. Mesh data in cartesian coordinates
and particle data can be read by this frontend.

To load the first in-file iteration of a openPMD datasets using the standard HDF5
output format:

.. code-block:: python

   import yt
   ds = yt.load('example-3d/hdf5/data00000100.h5')

If you operate on large files, you may want to modify the virtual chunking behaviour through
``open_pmd_virtual_gridsize``. The supplied value is an estimate of the size of a single read request
for each particle attribute/mesh (in Byte).

.. code-block:: python

  import yt
  ds = yt.load('example-3d/hdf5/data00000100.h5', open_pmd_virtual_gridsize=10e4)
  sp = yt.SlicePlot(ds, 'x', 'rho')
  sp.show()

Particle data is fully supported:

.. code-block:: python

  import yt
  ds = yt.load('example-3d/hdf5/data00000100.h5')
  ad = f.all_data()
  ppp = yt.ParticlePhasePlot(ad, 'particle_position_y', 'particle_momentum_y', 'particle_weighting')
  ppp.show()

.. rubric:: Caveats

* 1D, 2D and 3D data is compatible, but lower dimensional data might yield
  strange results since it gets padded and treated as 3D. Extraneous dimensions are
  set to be of length 1.0m and have a width of one cell.
* The frontend has hardcoded logic for renaming the openPMD ``position``
  of particles to ``positionCoarse``

.. _loading-pyne-data:

PyNE Data
---------

`PyNE <http://pyne.io/>`_ is an open source nuclear engineering toolkit
maintained by the PyNE development team (pyne-dev@googlegroups.com).
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

yt will attempt to guess the fields in the file. For more control over the hydro fields or the particle fields, see :ref:`loading-ramses-data-args`.

yt also support the new way particles are handled introduced after
version ``stable_17_09`` (the version introduced after the 2017 Ramses
User Meeting). In this case, the file ``part_file_descriptor.txt``
containing the different fields in the particle files will be read. If
you use a custom version of RAMSES, make sure this file is up-to-date
and reflects the true layout of the particles.

yt supports outputs made by the mainline ``RAMSES`` code as well as the
``RAMSES-RT`` fork. Files produces by ``RAMSES-RT`` are recognized as such
based on the presence of a ``info_rt_*.txt`` file in the output directory.

.. note::
   for backward compatibility, particles from the
   ``part_XXXXX.outYYYYY`` files have the particle type ``io`` by
   default (including dark matter, stars, tracer particles, ...). Sink
   particles have the particle type ``sink``.

.. _loading-ramses-data-args:

Arguments passed to the load function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is possible to provide extra arguments to the load function when loading RAMSES datasets. Here is a list of the ones specific to RAMSES:

``fields``
    A list of fields to read from the hydro files. For example, in a pure
    hydro simulation with an extra custom field named ``my-awesome-field``, one
    would specify the fields argument following this example:

      .. code-block:: python

          import yt
          fields = ["Density",
                    "x-velocity", "y-velocity", "z-velocity",
                    "Pressure", "my-awesome-field"]
          ds = yt.load('output_00123/info_00123.txt', fields=fields)
          'my-awesome-field' in ds.field_list  # is True


``extra_particle_fields``
      A list of tuples describing extra particles fields to read in. By
      default, yt will try to detect as many fields as possible,
      assuming the extra ones to be double precision floats. This
      argument is useful if you have extra fields besides the particle mass,
      position, and velocity fields that yt cannot detect automatically. For
      example, for a dataset containing two extra particle integer fields named
      ``family`` and ``info``, one would do:

      .. code-block:: python

          import yt
          extra_fields = [('family', 'I'), ('info', 'I')]
          ds = yt.load("output_00001/info_00001.txt", extra_particle_fields=extra_fields)
          # ('all', 'family') and ('all', 'info') now in ds.field_list

      The format of the ``extra_particle_fields`` argument is as follows:
      ``[('field_name_1', 'type_1'), ..., ('field_name_n', 'type_n')]`` where
      the second element of the tuple follows the `python struct format
      convention
      <https://docs.python.org/3.5/library/struct.html#format-characters>`_.
      Note that if ``extra_particle_fields`` is defined, yt will not assume
      that the ``particle_birth_time`` and ``particle_metallicity`` fields
      are present in the dataset. If these fields are present, they must be
      explicitly enumerated in the ``extra_particle_fields`` argument.

``cosmological``
      Force yt to consider a simulation to be cosmological or
      not. This may be useful for some specific simulations e.g. that
      run down to negative redshifts.

``bbox``
      The subbox to load. yt will only read CPUs intersecting with the
      subbox. This is especially useful for large simulations or
      zoom-in simulations, where you don't want to have access to data
      outside of a small region of interest. This argument will prevent
      yt from loading AMR files outside the subbox and will hence
      spare memory and time.
      For example, one could use

      .. code-block:: python

          import yt
          # Only load a small cube of size (0.1)**3
          bbox = [[0., 0., 0.], [0.1, 0.1, 0.1]]
          ds = yt.load('output_00001/info_00001.txt', bbox=bbox)

          # See the note below for the following examples
          ds.right_edge == [1, 1, 1]             # is True

          ad = ds.all_data()
          ad['particle_position_x'].max() > 0.1  # _may_ be True

          bb = ds.box(left_edge=bbox[0], right_edge=bbox[1])
          bb['particle_position_x'].max() < 0.1  # is True

      .. note::
         When using the bbox argument, yt will read all the CPUs
         intersecting with the subbox. However it may also read some
         data *outside* the selected region. This is due to the fact
         that domains have a complicated shape when using Hilbert
         ordering. Internally, yt will hence assume the loaded dataset
         covers the entire simulation. If you only want the data from
         the selected region, you may want to use ``ds.box(...)``.

      .. note::
         The ``bbox`` feature is only available for datasets using
         Hilbert ordering.

Adding custom particle fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three way to make yt detect all the particle fields. For example, if you wish to make yt detect the birth time and metallicity of your particles, use one of these methods

1. ``yt.load`` method. Whenever loading a dataset, add the extra particle fields as a keyword argument to the ``yt.load`` call.

   .. code-block:: python

      import yt
      epf = [('particle_birth_time', 'd'), ('particle_metallicity', 'd')]
      ds = yt.load('dataset', extra_particle_fields=epf)

      ('io', 'particle_birth_time') in ds.derived_field_list  # is True
      ('io', 'particle_metallicity') in ds.derived_field_list  # is True

2. yt config method. If you don't want to pass the arguments for each call of ``yt.load``, you can add in your configuration

   .. code-block:: none

      [ramses-particles]
      fields = particle_position_x, d
               particle_position_y, d
               particle_position_z, d
               particle_velocity_x, d
               particle_velocity_y, d
               particle_velocity_z, d
               particle_mass, d
               particle_identifier, i
               particle_refinement_level, I
               particle_birth_time, d
               particle_metallicity, d

3. New RAMSES way. Recent versions of RAMSES automatically write in their output an ``hydro_file_descriptor.txt`` file that gives information about which field is where. If you wish, you can simply create such a file in the folder containing the ``info_xxxxx.txt`` file

   .. code-block:: none

      # version:  1
      # ivar, variable_name, variable_type
       1, position_x, d
       2, position_y, d
       3, position_z, d
       4, velocity_x, d
       5, velocity_y, d
       6, velocity_z, d
       7, mass, d
       8, identity, i
       9, levelp, i
      10, birth_time, d
      11, metallicity, d

   It is important to note that this file should not end with an empty line (but in this case with ``11, metallicity, d``).

.. note::

   The kind (``i``, ``d``, ``I``, ...) of the field follow the `python convention <https://docs.python.org/3.5/library/struct.html#format-characters>`_.



Customizing the particle type association
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In verions of RAMSES more recent than December 2017, particles carry
along a ``family`` array. The value of this array gives the kind of
the particle, e.g. 1 for dark matter. It is possible to customize the
association between particle type and family by customizing the yt
config (see :ref:`configuration-file`), adding

.. code-block:: none

   [ramses-families]
   gas_tracer = 100
   star_tracer = 101
   dm = 0
   star = 1



Particle ages and formation times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For non-cosmological simulations, particle ages are stored in physical units on
disk. To access the birth time for the particles, use the
``particle_birth_time`` field. The time recorded in this field is relative to
the beginning of the simulation. Particles that were present in the initial
conditions will have negative values for ``particle_birth_time``.

For cosmological simulations that include star particles, RAMSES stores particle
formation times as conformal times. To access the formation time field data in
conformal units use the ``conformal_birth_time`` field. This will return the
formation times of particles in the simulation in conformal units as a
dimensionless array. To access the formation time in physical units, use the
``particle_birth_time`` field. Finally, to access the ages of star particles in
your simulation, use the ``star_age`` field. Note that this field is defined for
all particle types but will only make sense for star particles.

For simulations conducted in Newtownian coordinates, with no cosmology or
comoving expansion, the time is equal to zero at the beginning of the
simulation. That means that particles present in the initial conditions may have
negative birth times. This can happen, for example, in idealized isolated galaxy
simulations, where star particles are included in the initial conditions. For
simulations conducted in cosmological comoving units, the time is equal to zero
at the big bang, and all particles should have positive values for the
``particle_birth_time`` field.

To help clarify the above discussion, the following table describes the meaning
of the various particle formation time and age fields:

+------------------+--------------------------+--------------------------------+
| Simulation type  | Field name               | Description                    |
|==================+==========================+================================+
| cosmological     | ``conformal_birth_time`` | Formation time in conformal    |
|                  |                          | units (dimensionless)          |
+------------------+--------------------------+--------------------------------+
| any              | ``particle_birth_time``  | The time relative to the       |
|                  |                          | beginning of the simulation    |
|                  |                          | when the particle was formed.  |
|                  |                          | For non-cosmological           |
|                  |                          | simulations, this field will   |
|                  |                          | have positive values for       |
|                  |                          | particles formed during the    |
|                  |                          | simulation and negative for    |
|                  |                          | particles of finite age in the |
|                  |                          | initial conditions. For        |
|                  |                          | cosmological simulations this  |
|                  |                          | is the time the particle       |
|                  |                          | formed relative to the big     |
|                  |                          | bang, therefore the value of   |
|                  |                          | this field should be between   |
|                  |                          | 0 and 13.7 Gyr.                |
+------------------+--------------------------+--------------------------------+
| any              | ``star_age``             | Age of the particle.           |
|                  |                          | Only physically meaningful for |
|                  |                          | stars and particles that       |
|                  |                          | formed dynamically during the  |
|                  |                          | simulation.                    |
+------------------+--------------------------+--------------------------------+

RAMSES datasets produced by a version of the code newer than November 2017
contain the metadata necessary for yt to automatically distinguish between star
particles and other particle types. If you are working with a dataset produced
by a version of RAMSES older than November 2017, yt will only automatically
recognize a single particle ``io``. It may be convenient to define a particle
filter in your scripts to distinguish between particles present in the initial
conditions and particles that formed dynamically during the simulation by
filtering particles with ``"conformal_birth_time"`` values equal to zero and not
equal to zero.  An example particle filter definition for dynamically formed
stars might look like this:

.. code-block:: python

    @yt.particle_filter(requires=["conformal_birth_time"],
                        filtered_type='io')
    def stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "conformal_birth_time"] != 0
        return filter

For a cosmological simulation, this filter will distinguish between stars and
dark matter particles.
        
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

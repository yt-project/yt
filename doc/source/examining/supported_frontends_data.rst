.. _loading-data-from-supported-codes:

Loading Data from Supported Codes
=================================

This section contains information on how to load data into ``yt`` from
supported codes, as well as some important caveats about different
data formats.

.. _loading-enzo-data:

Enzo Data
---------

Enzo data is fully supported and cared for by Matthew Turk.  To load an Enzo
dataset, you can use the ``load`` command provided by ``yt.mods`` and supply to
it the parameter file name.  This would be the name of the output file, and it
contains no extension.  For instance, if you have the following files:

.. code-block:: none

   DD0010/
   DD0010/data0010
   DD0010/data0010.hierarchy
   DD0010/data0010.cpu0000
   DD0010/data0010.cpu0001
   DD0010/data0010.cpu0002
   DD0010/data0010.cpu0003

You would feed the ``load`` command the filename ``DD0010/data0010`` as
mentioned.

.. code-block:: python

   from yt.mods import *
   pf = load("DD0010/data0010")

.. rubric:: Caveats

* There are no major caveats for Enzo usage
* Units should be correct, if you utilize standard unit-setting routines.  yt
  will notify you if it cannot determine the units, although this
  notification will be passive.
* 2D and 1D data are supported, but the extraneous dimensions are set to be
  of length 1.0

.. _loading-orion-data:

Orion Data
----------

Orion data is fully supported. To load an Orion dataset, you can use the
``load`` command provided by ``yt.mods`` and supply to it the directory file
name.  **You must also have the ``inputs`` file in the base directory.** For
instance, if you were in a directory with the following files:

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

   from yt.mods import *
   pf = load("pltgmlcs5600")

.. rubric:: Caveats

* There are no major caveats for Orion usage
* Star particles are not supported at the current time

.. _loading-flash-data:

FLASH Data
----------

FLASH HDF5 data is fully supported and cared for by John ZuHone.  To load a
FLASH dataset, you can use the ``load`` command provided by ``yt.mods`` and
supply to it the file name of a plot file or checkpoint file.  Particle
files are not currently directly loadable by themselves, due to the
fact that they typically lack grid information. For instance, if you were in a directory with
the following files:

.. code-block:: none

   cosmoSim_coolhdf5_chk_0026

You would feed it the filename ``cosmoSim_coolhdf5_chk_0026``:

.. code-block:: python

   from yt.mods import *
   pf = load("cosmoSim_coolhdf5_chk_0026")

If you have a FLASH particle file that was created at the same time as
a plotfile or checkpoint file (therefore having particle data
consistent with the grid structure of the latter), its data may be loaded with the
``particle_filename`` optional argument:

.. code-block:: python

    from yt.mods import *
    pf = load("radio_halo_1kpc_hdf5_plt_cnt_0100", particle_filename="radio_halo_1kpc_hdf5_part_0100")

.. rubric:: Caveats

* Please be careful that the units are correctly utilized; yt assumes cgs
* Velocities and length units will be scaled to comoving coordinates if yt is
  able to discern you are examining a cosmology simulation; particle and grid
  positions will not be.
* Domains may be visualized assuming periodicity.

Athena Data
-----------

Athena 4.x VTK data is *mostly* supported and cared for by John
ZuHone. Both uniform grid and SMR datasets are supported. 

Loading Athena datasets is slightly different depending on whether
your dataset came from a serial or a parallel run. If the data came
from a serial run or you have joined the VTK files together using the
Athena tool ``join_vtk``, you can load the data like this:

.. code-block:: python

   from yt.mods import *
   pf = load("kh.0010.vtk")

The filename corresponds to the file on SMR level 0, whereas if there
are multiple levels the corresponding files will be picked up
automatically, assuming they are laid out in ``lev*`` subdirectories
under the directory where the base file is located.

For parallel datasets, yt assumes that they are laid out in
directories named ``id*``, one for each processor number, each with
``lev*`` subdirectories for additional refinement levels. To load this
data, call ``load`` with the base file in the ``id0`` directory:

.. code-block:: python

   from yt.mods import *
   pf = load("id0/kh.0010.vtk")

which will pick up all of the files in the different ``id*`` directories for
the entire dataset. 

yt works in cgs ("Gaussian") units, but Athena data is not
normally stored in these units. If you would like to convert data to
cgs units, you may supply conversions for length, time, and density to ``load``:

.. code-block:: python

   from yt.mods import *
   pf = load("id0/cluster_merger.0250.vtk", 
          parameters={"LengthUnits":3.0856e24,
                               "TimeUnits":3.1557e13,"DensityUnits":1.67e-24)

This means that the yt fields (e.g. ``Density``, ``x-velocity``,
``Bx``) will be in cgs units, but the Athena fields (e.g.,
``density``, ``velocity_x``, ``cell_centered_B_x``) will be in code
units. 

.. rubric:: Caveats

* yt primarily works with primitive variables. If the Athena
  dataset contains conservative variables, the yt primitive fields will be generated from the
  conserved variables on disk. 
* Domains may be visualized assuming periodicity.
* Particle list data is currently unsupported.
* In some parallel Athena datasets, it is possible for a grid from one
  refinement level to overlap with more than one grid on the parent
  level. This may result in unpredictable behavior for some analysis
  or visualization tasks. 

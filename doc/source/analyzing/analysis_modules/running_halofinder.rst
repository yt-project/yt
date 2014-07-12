.. _halo_finding:

Halo Finding
============
.. sectionauthor:: Stephen Skory <sskory@physics.ucsd.edu>

There are four methods of finding particle haloes in yt. The recommended and default method is called HOP, a 
method described in `Eisenstein and Hut (1998) <http://adsabs.harvard.edu/abs/1998ApJ...498..137E>`_. 
A basic friends-of-friends (e.g. `Efstathiou et al. (1985) <http://adsabs.harvard.edu/abs/1985ApJS...57..241E>`_)
halo finder is also implemented.
Parallel HOP (`Skory et al. (2010) <http://adsabs.harvard.edu/abs/2010ApJS..191...43S>`_)
is a true parallelization of the HOP method can analyze massive datasets on
hundreds of processors.
Finally Rockstar (`Behroozi et a. (2011) <http://adsabs.harvard.edu/abs/2011arXiv1110.4372B>`_)
is a 6D-phase space halo finder developed by Peter Behroozi
that excels in finding subhalos and substrcture,
but does not allow multiple particle masses.

HOP
---

The version of HOP used in yt is an upgraded version of the `publicly available HOP code 
<http://cmb.as.arizona.edu/~eisenste/hop/hop.html>`_. Support for 64-bit floats and integers has been
added, as well as parallel analysis through spatial decomposition. HOP builds groups in this fashion:

  1. Estimates the local density at each particle using a smoothing kernel.
  2. Builds chains of linked particles by 'hopping' from one particle to its densest neighbor.
     A particle which is its own densest neighbor is the end of the chain.
  3. All chains that share the same densest particle are grouped together.
  4. Groups are included, linked together, or discarded depending on the user-supplied over density
     threshold parameter. The default is 160.0.

Please see the `HOP method paper <http://adsabs.harvard.edu/abs/1998ApJ...498..137E>`_ 
for full details.

Friends-of-Friends
------------------

The version of FoF in yt is based on the `publicly available FoF code <http://www-hpcc.astro.washington.edu/tools/fof.html>`_ from the University of Washington. Like HOP,
FoF supports parallel analysis through spatial decomposition. FoF is much simpler than HOP:

  1. From the total number of particles, and the volume of the region, the average
     inter-particle spacing is calculated.
  2. Pairs of particles closer together than some fraction of the average inter-particle spacing
     (the default is 0.2) are linked together. Particles can be paired with more than one other particle.
  3. The final groups are formed the networks of particles linked together by friends, hence the name.

.. warning:: The FoF halo finder in yt is not thoroughly tested! It is probably fine to use, but you
   are strongly encouraged to check your results against the data for errors.

Running HaloFinder
------------------

Running HOP on a dataset is straightforward

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  halo_list = HaloFinder(pf)

Running FoF is similar:

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  halo_list = FOFHaloFinder(pf)

Halo Data Access
----------------

``halo_list`` is a list of ``Halo`` class objects ordered by decreasing halo mass. A ``Halo`` object
has convenient ways to access halo data. This loop will print the location of the center of mass
for each halo found

.. code-block:: python

  for halo in halo_list:
      print halo.center_of_mass()

All the methods are:

  * .center_of_mass() - the center of mass for the halo.
  * .maximum_density() - the maximum density in "HOP" units.
  * .maximum_density_location() - the location of the maximum density particle in the HOP halo.
  * .total_mass() - the mass of the halo in Msol (not Msol/h).
  * .bulk_velocity() - the velocity of the center of mass of the halo in simulation units.
  * .maximum_radius() - the distance from the center of mass to the most distant particle in the halo
    in simulation units.
  * .get_size() - the number of particles in the halo.
  * .get_sphere() - returns an an EnzoSphere object using the center of mass and maximum radius.
  * .virial_mass(virial_overdensity=float, bins=int) - Finds the virial
    mass for a halo using just the particles. This is inferior to the full
    Halo Profiler extension (:ref:`halo_profiling`), but useful nonetheless in some cases.
    Returns the mass in Msol, or -1 if the halo is not virialized.
    Defaults: ``virial_overdensity=200.0`` and ``bins=300``.
  * .virial_radius(virial_overdensity=float, bins=int) - Fins the virial
    radius of the halo using just the particles. Returns the radius in code
    units, or -1 if the halo is not virialized.
    Defaults: ``virial_overdensity=200.0`` and ``bins=300``.

.. note:: For FOF the maximum density value is meaningless and is set to -1 by default. For FOF
   the maximum density location will be identical to the center of mass location.

For each halo the data for the particles in the halo can be accessed like this

.. code-block:: python

  for halo in halo_list:
      print halo["particle_index"]
      print halo["particle_position_x"] # in simulation units

Halo List Data Access
---------------------

These are methods that operate on the list of halo objects, rather than on the
haloes themselves (e.g. ``halo_list.write_out()`` instead of ``halo_list[0].center_of_mass()``).
For example, The command

.. code-block:: python

  halo_list.write_out("HaloAnalysis.out")

will output the haloes to a text file named ``HaloAnalysis.out``.

  * .write_out(``name``) - Writes out the center of mass, maximum density point,
    number of particles, mass, index, bulk velocity and maximum radius for all the haloes
    to a text file ``name``.
  * .write_particle_lists(``name``) - Writes the data for the particles in haloes
    (position, velocity, mass and particle index) to a HDF5 file with prefix ``name``, or one HDF5
    file per CPU when running in parallel.
  * .write_particle_lists_txt(``name``) - Writes out one text file with prefix ``name`` that gives the
    location of the particle data for haloes in the HDF5 files. This is only
    necessary when running in parallel.
  * .dump(``basename``) - Calls all of the above three functions using 
    ``basename`` in each. This function is meant to be used in combination with
    loading halos off disk (:ref:`load_haloes`).
  * .nearest_neighbors_3D(haloID, num_neighbors=int, search_radius=float) - 
    For a given halo ``haloID``, this finds the ``num_neighbors`` nearest (periodic)
    neighbors that are within ``search_radius`` distance from it.
    It returns a list of the neighbors distances and ID with format
    [distance,haloID]. Defaults: ``num_neighbors=7``, ``search_radius=0.2``.
  * .nearest_neighbors_2D(haloID, num_neighbors=int, search_radius=float, proj_dim={0,1,2}) -
    Similarly to the 3D search, this finds the nearest (periodic) neighbors to a halo, but
    with the positions of the haloes projected onto a 2D plane. The normal to the
    projection plane is set with ``proj_dim``, which is set to {0,1,2} for the
    {x,y,z}-axis. Defaults: ``num_neighbors=7``, ``search_radius=0.2`` and ``proj_dim=0``.
    Returns a list of neighbors in the same format as the 3D case, but the distances
    are the 2D projected distance.

.. _load_haloes:

Loading Haloes Off Disk
-----------------------

It is possible to load haloes off disk and use them as if they had just been
located by the halo finder. This has at least two advantages.  Quite obviously
this means that if the halos are properly saved (e.g. ``haloes.dump()``, see
above and below), halo finding does not need to be run again, saving time.
Another benefit is loaded haloes only use as much memory as needed because the
particle data for the haloes is loaded off disk on demand. If only a few haloes
are being examined, a dataset that required parallel analysis for halo finding
can be analyzed in serial, interactively.

The first step is to save the haloes in a consistent manner, which is made
simple with the ``.dump()`` function:

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  haloes = HaloFinder(pf)
  haloes.dump("basename")

It is easy to load the halos using the ``LoadHaloes`` class:

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  haloes = LoadHaloes(pf, "basename")

Everything that can be done with ``haloes`` in the first example should be
possible with ``haloes`` in the second.

General Parallel Halo Analysis
------------------------------

Both the HOP and FoF halo finders can run in parallel using simple spatial decomposition.
In order to run them
in parallel it is helpful to understand how it works.

Below in the first plot (i) is a simplified depiction of three haloes labeled 1,2 and 3:

.. image:: _images/ParallelHaloFinder.png
   :width: 500

Halo 3 is twice reflected around the periodic boundary conditions.

In (ii), the volume has been
sub-divided into four equal subregions, A,B,C and D, shown with dotted lines. Notice that halo 2
is now in two different subregions,
C and D, and that halo 3 is now in three, A, B and D. If the halo finder is run on these four separate subregions,
halo 1 is be identified as a single halo, but haloes 2 and 3 are split up into multiple haloes, which is incorrect.
The solution is to give each subregion padding to oversample into neighboring regions.

In (iii), subregion C has oversampled into the other three regions, with the periodic boundary conditions taken
into account, shown by dot-dashed lines. The other subregions oversample in a similar way.

The halo finder is then run on each padded subregion independently and simultaneously.
By oversampling like this, haloes 2 and 3 will both be enclosed fully in at least one subregion and
identified completely.

Haloes identified with centers of mass inside the padded part of a subregion are thrown out, eliminating
the problem of halo duplication. The centers for the three haloes are shown with stars. Halo 1 will
belong to subregion A, 2 to C and 3 to B.

Parallel HaloFinder padding
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run with parallel halo finding, there is a slight modification to the script

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  halo_list = HaloFinder(pf,padding=0.02)
  # --or--
  halo_list = FOFHaloFinder(pf,padding=0.02)

The ``padding`` parameter is in simulation units and defaults to 0.02. This parameter is how much padding
is added to each of the six sides of a subregion. This value should be 2x-3x larger than the largest
expected halo in the simulation. It is unlikely, of course, that the largest object in the simulation
will be on a subregion boundary, but there is no way of knowing before the halo finder is run.

In general, a little bit of padding goes a long way, and too much just slows down the analysis and doesn't
improve the answer (but doesn't change it). 
It may be worth your time to run the parallel halo finder at a few paddings to
find the right amount, especially if you're analyzing many similar datasets.

Parallel HOP
------------

**Parallel HOP** (not to be confused with HOP running in parallel as described
above) is a wholly-new halo finder based on the HOP method.
For extensive details and benchmarks of Parallel HOP, please see the
pre-print version of the `method paper <http://adsabs.harvard.edu/abs/2010ApJS..191...43S>`_ at
arXiv.org.
While the method
of parallelization described above can be quite effective, it has its limits.
In particular
for highly unbalanced datasets, where most of the particles are in a single
part of the simulation's volume, it can become impossible to subdivide the
volume sufficiently to fit a subvolume into a single node's memory.

Parallel HOP is designed to be parallel at all levels of operation. There is
a minimal amount of copied data across tasks. Unlike the parallel method above,
whole haloes do not need to exist entirely in a single subvolume. In fact, a
halo may have particles in several subvolumes simultaneously without a problem.

Parallel HOP is appropriate for very large datasets where normal HOP, or
the parallel method described above, won't work. For smaller datasets, it is
actually faster to use the simpler methods above because the mechanisms employed for
full parallelism are somewhat expensive.
Whether to use Parallel HOP or not depends on the number of particles and
the size of the largest object in the simulation.
Because the padding of the other parallel method described above depends on
the relative size to the box of the largest object, for smaller cosmologies
that method may not work.
If the largest object is quite large, the minimum padding will be a
significant fraction of the full volume, and therefore the minimum number of
particles per task can stay quite high.
Below and including 256^3 particles, the other parallel methods are likely
faster.
However, above this and for smaller cosmologies (100 Mpc/h and smaller),
Parallel HOP will offer better performance.

The haloes identified by Parallel HOP are slightly different than normal HOP
when run on the same dataset with the same over-density threshold.
For a given threshold value, a few haloes have slightly different numbers of particles.
Overall, it is not a big difference. In fact, changing the threshold value by
a percent gives a far greater difference than the differences between HOP and
Parallel HOP.

HOP and Parallel HOP both use `KD Trees <http://en.wikipedia.org/wiki/Kd_tree>`_
for nearest-neighbor searches.
Parallel HOP uses the Fortran version of
`KDTREE 2 <http://arxiv.org/abs/physics/0408067>`_ written by Matthew B. Kennel.
The KD Tree in normal HOP calculates the distances
between particles incorrectly by approximately one part in a million.
KDTREE 2 is far more accurate (up to machine error),
and this slight difference is sufficient to make perfect agreement between
normal and Parallel HOP impossible.
Therefore Parallel HOP is not a direct substitution for
normal HOP, but is very similar.

Running Parallel HOP
^^^^^^^^^^^^^^^^^^^^

Note: This is probably broken now that the Fortran kdtree has been removed.

In the simplest form, Parallel HOP is run very similarly to the other halo finders.
In the example below, Parallel HOP will be run on a dataset with all the default
values. Parallel HOP can be run in serial, but as mentioned above, it is
slower than normal HOP.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  halo_list = parallelHF(pf)

Parallel HOP has these user-set options:

  * ``threshold``, positive float: This is the same as the option for normal HOP. Default=160.0.
  * ``dm_only``, True/False: Whether or not to include particles other than dark
    matter when building haloes. Default=True.
  * ``resize``, True/False: Parallel HOP can load-balance the particles, such that
    each subvolume has the same number of particles.
    In general, this option is a good idea for simulations' volumes
    smaller than about 300 Mpc/h, and absolutely required for those under
    100 Mpc/h. For larger volumes the particles are distributed evenly enough
    that this option is unnecessary. Default=True.
  * ``sample``, positive float: In order to load-balance, a random subset of the
    particle positions are read off disk, and the load-balancing routine is
    applied to them. This parameter controls what fraction of the full dataset
    population is used. Larger values result in more accurate load-balancing,
    and smaller values are faster. The value cannot be too large as the data
    for the subset of particles is communicated to one task for
    load-balancing (meaning a value
    of 1.0 will not work on very large datasets).
    Tests show that values as low as 0.0003 keep the min/max variation between
    tasks below 10%. Default = 0.03.
  * ``rearrange``, True/False: The KD Tree used by Parallel HOP can make an
    internal copy of the particle data which increases the speed of nearest
    neighbor searches by approximately 20%. The only reason to turn this option
    off is if memory is a concern. Default=True.
  * ``safety``, positive float: Unlike the simpler parallel method, Parallel
    HOP calculates the padding automatically. The padding is a
    function of the inter-particle spacing inside each subvolume. This parameter
    is multiplied to the padding distance to increase the padding volume to account for
    density variations on the boundaries of the subvolumes. Increasing this
    parameter beyond a certain point will have no effect other than consuming
    more memory and slowing down runtimes.
    Reducing it will speed up the calculation and use less memory, but 
    going too far will result in degraded halo finding.
    Default=1.5, but values as low as 1.0 will probably work for many datasets.
  * ``fancy_padding``, True/False: When this is set to True, the amount of padding
    is calculated independently for each of the six faces of each subvolume. When this is
    False, the padding is the same on all six faces. There is generally no
    good reason to set this to False. Default=True.
  * ``premerge``, True/False: This option will pre-merge only the most dense
    haloes in each subvolume, before haloes are merged on the global level. In
    some cases this can speed up the runtime by a factor of two and reduce peak memory
    greatly. At worst it slows down the runtime by a small amount. It has the
    side-effect of changing the haloes slightly as a function of task count. Put in
    another way, two otherwise identical runs of Parallel HOP on a dataset will end
    up with very slightly different haloes when run with two different task counts
    with this option turned on.  Not all haloes are changed between runs.  This is
    due to the way merging happens in HOP - pre-merging destroys the global
    determinacy of halo merging. Default=True.
  * ``tree``, string: There are two kD-trees that may be used as part of the
    halo-finding process. The Fortran ("F") one is (presently) faster, but requires
    more memory. One based on `scipy.spatial
    <http://docs.scipy.org/doc/scipy/reference/spatial.html>`_ utilizes
    Cython ("C") and is (presently) slower, but is more memory efficient.
    Default = "F".

All the same halo data can be accessed from Parallel HOP haloes as with the other halo finders.
However, when running in parallel, there are some
important differences in the output of a couple of these functions.

  * .write_particle_lists(``name``) - Because haloes may exist in more than
    one subvolume, particle data for a halo may be saved in more than one HDF5 file.
  * .write_particle_lists_txt(``name``) - If the particles for a halo is saved
    in more than one HDF5 file, there will be more than one HDF5 file listed for
    each halo in the text file.

In this example script below, Parallel HOP is run on a dataset and the results
saved to files. The summary of the haloes to ``ParallelHopAnalysis.out``, the
particles to files named ``parts????.h5`` and the list of haloes in HDF5 files
to ``parts.txt``.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load("data0001")
  halo_list = parallelHF(pf, threshold=80.0, dm_only=True, resize=False, 
  rearrange=True, safety=1.5, premerge=True)
  halo_list.write_out("ParallelHopAnalysis.out")
  halo_list.write_particle_list("parts")
  halo_list.write_particle_lists_txt("parts")

Halo Finding In A Subvolume
---------------------------

It is possible to run any of the halo finders over a subvolume.
This may be advantageous when only one object or region of a simulation
is being analyzed.
The subvolume must be a ``region`` and cannot be a
non-rectilinear shape.
The halo finding can be performed in parallel on a subvolume, but it may
not be necessary depending on the size of the subvolume.
Below is a simple example for HOP; the other halo finders use the same
``subvolume`` keyword identically.

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.api import *
  pf = load('data0458')
  # Note that the first term below, [0.5]*3, defines the center of
  # the region and is not used. It can be any value.
  sv = pf.region([0.5]*3, [0.21, .21, .72], [.28, .28, .79])
  halos = HaloFinder(pf, subvolume = sv)
  halos.write_out("sv.out")


Rockstar Halo Finding
=====================
.. sectionauthor:: Matthew Turk <matthewturk@gmail.com>
.. sectionauthor:: Christopher Erick Moody<cemoody@ucsc.edu>
.. sectionauthor:: Stephen Skory <s@skory.us>

Rockstar uses an adaptive hierarchical refinement of friends-of-friends 
groups in six phase-space dimensions and one time dimension, which 
allows for robust (grid-independent, shape-independent, and noise-
resilient) tracking of substructure. The code is prepackaged with yt, 
but also `separately available <http://code.google.com/p/rockstar>`_. The lead 
developer is Peter Behroozi, and the methods are described in `Behroozi
et al. 2011 <http://rockstar.googlecode.com/files/rockstar_ap101911.pdf>`_. 

.. note:: At the moment, Rockstar does not support multiple particle masses, 
  instead using a fixed particle mass. This will not affect most dark matter 
  simulations, but does make it less useful for finding halos from the stellar
  mass. Also note that halo finding in a subvolume is not supported by
  Rockstar.

To run the Rockstar Halo finding, you must launch python with MPI and 
parallelization enabled. While Rockstar itself does not require MPI to run, 
the MPI libraries allow yt to distribute particle information across multiple 
nodes.

.. warning:: At the moment, running Rockstar inside of yt on multiple compute nodes
   connected by an Infiniband network can be problematic. Therefore, for now
   we recommend forcing the use of the non-Infiniband network (e.g. Ethernet)
   using this flag: ``--mca btl ^openib``.
   For example, here is how Rockstar might be called using 24 cores:
   ``mpirun -n 24 --mca btl ^openib python ./run_rockstar.py --parallel``.

Designing the python script itself is straightforward:

.. code-block:: python

  from yt.mods import *
  from yt.analysis_modules.halo_finding.rockstar.api import RockstarHaloFinder

  #find all of our simulation files
  files = glob.glob("Enzo_64/DD*/\*index")
  #hopefully the file name order is chronological
  files.sort()
  ts = DatasetSeries.from_filenames(files[:])
  rh = RockstarHaloFinder(ts)
  rh.run()

The script above configures the Halo finder, launches a server process which 
disseminates run information and coordinates writer-reader processes. 
Afterwards, it launches reader and writer tasks, filling the available MPI 
slots, which alternately read particle information and analyze for halo 
content.

The RockstarHaloFinder class has these options:
  * ``dm_type``, the index of the dark matter particle. Default is 1. 
  * ``outbase``, This is where the out*list files that Rockstar makes should be
    placed. Default is 'rockstar_halos'.
  * ``num_readers``, the number of reader tasks (which are idle most of the 
    time.) Default is 1.
  * ``num_writers``, the number of writer tasks (which are fed particles and
    do most of the analysis). Default is MPI_TASKS-num_readers-1. 
    If left undefined, the above options are automatically 
    configured from the number of available MPI tasks.
  * ``force_res``, the resolution that Rockstar uses for various calculations
    and smoothing lengths. This is in units of Mpc/h.
    If no value is provided, this parameter is automatically set to
    the width of the smallest grid element in the simulation from the
    last data snapshot (i.e. the one where time has evolved the
    longest) in the time series:
    ``pf_last.index.get_smallest_dx() * pf_last['mpch']``.
  * ``total_particles``, if supplied, this is a pre-calculated
    total number of dark matter
    particles present in the simulation. For example, this is useful
    when analyzing a series of snapshots where the number of dark
    matter particles should not change and this will save some disk
    access time. If left unspecified, it will
    be calculated automatically. Default: ``None``.
  * ``dm_only``, if set to ``True``, it will be assumed that there are
    only dark matter particles present in the simulation.
    This option does not modify the halos found by Rockstar, however
    this option can save disk access time if there are no star particles
    (or other non-dark matter particles) in the simulation. Default: ``False``.


Output Analysis
---------------

Rockstar dumps halo information in a series of text (halo*list and 
out*list) and binary (halo*bin) files inside the ``outbase`` directory. 
We use the halo list classes to recover the information. 

Inside the ``outbase`` directory there is a text file named ``pfs.txt``
that records the connection between pf names and the Rockstar file names.

The halo list can be automatically generated from the RockstarHaloFinder 
object by calling ``RockstarHaloFinder.halo_list()``. Alternatively, the halo
lists can be built from the RockstarHaloList class directly 
``LoadRockstarHalos(pf,'outbase/out_0.list')``.

.. code-block:: python
    
    rh = RockstarHaloFinder(pf)
    #First method of creating the halo lists:
    halo_list = rh.halo_list()    
    #Alternate method of creating halo_list:
    halo_list = LoadRockstarHalos(pf, 'rockstar_halos/out_0.list')

The above ``halo_list`` is very similar to any other list of halos loaded off
disk.
It is possible to access particle data and use the halos in a manner like any
other halo object, and the particle data is only loaded on demand.
Additionally, each halo object has additional information attached that is
pulled directly from the Rockstar output:

.. code-block:: python

    >>> halo_list[0].supp
    Out[3]: 
    {'J': array([ -6.15271728e+15,  -1.36593609e+17,  -7.80776865e+16], dtype=float32),
     'bulkvel': array([-132.05046082,   11.53190422,   42.16183472], dtype=float32),
     'child_r': 2.6411054,
     'corevel': array([-132.05046082,   11.53190422,   42.16183472], dtype=float32),
     'desc': 0,
     'energy': -8.106986e+21,
     'flags': 1,
     'id': 166,
     'm': 1.5341227e+15,
     'mgrav': 1.5341227e+15,
     'min_bulkvel_err': 1821.8152,
     'min_pos_err': 0.00049575343,
     'min_vel_err': 1821.8152,
     'n_core': 1958,
     'num_child_particles': 2764,
     'num_p': 2409,
     'p_start': 6540,
     'pos': array([   0.20197368,    0.54656458,    0.11256824, -104.33285522,
             29.02485085,   43.5154953 ], dtype=float32),
     'r': 0.018403014,
     'rs': 0.0026318002,
     'rvmax': 1133.2,
     'spin': 0.035755754,
     'vmax': 1877.125,
     'vrms': 1886.2648}

Installation
------------

The Rockstar is slightly patched and modified to run as a library inside of 
yt. By default it will be built with yt using the ``install_script.sh``.
If it wasn't installed, please make sure that the installation setting
``INST_ROCKSTAR=1`` is defined in the ``install_script.sh`` and re-run
the installation script.

Rockstar Inline with Enzo
-------------------------

It is possible to run Rockstar inline with Enzo. Setting up
Enzo with inline yt is covered
`here <http://enzo-project.org/doc/user_guide/EmbeddedPython.html>`_.
It is not necessary to run Enzo with load balancing off to use Rockstar.
Here is an example ``user_script.py``:

.. code-block:: python

    from yt.mods import *
    from yt.analysis_modules.halo_finding.api import *
    from yt.config import ytcfg
    from yt.analysis_modules.halo_finding.rockstar.api import *
    
    def main():
        import enzo
        pf = EnzoDatasetInMemory()
        mine = ytcfg.getint('yt','__topcomm_parallel_rank')
        size = ytcfg.getint('yt','__topcomm_parallel_size')

        # Call rockstar.
        ts = DatasetSeries([pf])
        outbase = "./rockstar_halos_%04d" % pf['NumberOfPythonTopGridCalls']
        rh = RockstarHaloFinder(ts, num_readers = size,
            outbase = outbase)
        rh.run()
    
        # Load the halos off disk.
        fname = outbase + "/out_0.list"
        rhalos = LoadRockstarHalos(pf, fname)


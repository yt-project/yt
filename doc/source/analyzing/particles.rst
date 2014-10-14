Examining and Manipulating Particles
====================================

``yt`` has support for reading and manipulating particles.  You can access the
particles as you would any other data field; additionally, derived fields that
operate on particles can be added as would any other derived field, as long as
the parameter *particle_type* is set to ``True`` in the call to
:func:`add_field`.  However, with that, there are a few caveats.  Particle
support in ``yt`` is not by any means an afterthought, but it was developed
relatively late in comparison to baryon and field-based analysis, and is not as
mature.

.. note:: If you are having trouble with particles, email the mailing list!

Using Particles
---------------

Many particle operations can be conducted indirectly, which will serve to
reduce memory usage as well as handle any problems that might arise from
spatial selection of particles.

For instance, :class:`~yt.analysis_modules.halo_finding.halo_objects.Halo` 
objects have a number of operations that
can transparently calculate center of mass of particles, bulk velocity, and so
on.  Use those instead of obtaining the fields directly.  Furthermore, any of
the spatially-addressable objects described in :ref:`using-objects` will
automatically select particles based on the region of space they describe, and
the quantities (:ref:`derived-quantities`) in those objects will operate on
particle fields.

(For information on halo finding, see :ref:`cookbook-halo_finding`)

.. warning:: If you use the built-in methods of interacting with particles, you
             should be well off.  Otherwise, there are caveats!

Selection By Type
-----------------

Unfortunately, Enzo's mechanism for storing particle type is inconsistent.  The
parameter ``ParticleTypeInFile`` controls whether or not the field
``particle_type`` is written to disk; if it is set to 1, the field will be
written, but the default is 0 where the field is not written.  Without the
field ``particle_type`` the discriminator between particle types is exclusively
based on the field ``creation_time``.  Particles with ``creation_time`` greater
than 0.0 are star particles and those with ``creation_time`` equal to zero are
dark matter particles.

For simulations only including dark matter particles, this is not important, as
all of the particles will be of the same type.  However, selection of -- for
instance -- star particles in other simulations will require some care, and you
will need to do it differently depending on the value of
``ParticleTypeInFile``.

.. _selecting-creation-time:

Selecting Particles By Creation Time
++++++++++++++++++++++++++++++++++++

To select particles based on creation time, you must first create an index
array.  Python (and NumPy) allow indexing based on boolean values, so we will
do that.  Here is an example of selecting all star particles in the domain.

.. code-block:: python

   from yt.mods import *
   ds = load("galaxy1200.dir/galaxy1200")
   dd = ds.all_data()

   star_particles = dd["creation_time"] > 0.0
   print dd["ParticleMassMsun"][star_particles].max()
   print dd["ParticleMassMsun"][star_particles].min()
   print "Number of star particles", star_particles.sum()

Selecting Particles By Particle Type
++++++++++++++++++++++++++++++++++++

In Enzo, star particles are type 2.  So we will select using the boolean array
(as in :ref:`selecting-creation-time`) to select only the star particles.

.. code-block:: python

   from yt.mods import *
   ds = load("galaxy1200.dir/galaxy1200")
   dd = ds.all_data()

   star_particles = dd["particle_type"] == 2
   print dd["ParticleMassMsun"][star_particles].max()
   print dd["ParticleMassMsun"][star_particles].min()
   print "Number of star particles", star_particles.sum()

Memory
------

Unfortunately, as of right now, particle loading via spatially-selected objects
can be memory intensive.  The process that ``yt`` goes through to load
particles into memory in a 3D data object is to separate the grids into two
classes:

 * Fully-contained grids
 * Partially-contained grids

For the grids in the former category, the full set of particles residing in
those grids are loaded.  The ones in the second require that a
:class:`~yt.data_objects.data_containers.FakeGridForParticles` be created so 
that the particles residing in the region (as determined by their values of
``particle_position_x``, ``particle_position_y`` and ``particle_position_z``,
which must be loaded from disk) can be selected and cut from the full set of
particles.  This requires that the full position information for the particles
be loaded, which increases overall memory usage.

The Future
----------

The next version of ``yt`` will have a completely rewritten particle
infrastructure.  This version is currently in the testing phase, but has shown
to reduce memory overhead substantially as well as increase speed by a factor
of a few.  Both spatial selection (selection within an object) and selection by
type are extremely promising.

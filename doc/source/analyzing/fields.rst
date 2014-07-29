Fields in yt
============

The fundamental way to query data in yt is to access a field, either in its raw
form (by examining a data container) or a processed form (derived quantities,
projections, and so on).  "Field" is something of a loaded word, as it can
refer to quantities that are defined everywhere, which we refer to as "mesh" or
"fluid" fields, or discrete points that populate the domain, traditionally
thought of as "particle" fields.  The word "particle" here is gradually falling
out of favor, as these discrete fields can be any type of sparsely populated
data.

In previous versions of yt, there was a single mechanism of accessing fields on
a data container -- by their name, which was mandated to be a single string,
and which often varied between different code frontends.  yt 3.0 allows
for datasets containing multiple different types of fluid fields, mesh fields,
particles (with overlapping or disjoint lists of fields).  To enable accessing
these fields in a meaningful, simple way, the mechanism for accessing them has
changed to take an optional *field type* in addition to the *field name*.

As an example, we may be in a situation where have multiple types of particles
which possess the ``particle_position`` field.  In the case where a data
container, here called ``ad`` (short for "all data") contains a field, we can
specify which particular particle type we want to query:

.. code-block:: python

   print ad["humans", "particle_position"]
   print ad["dogs", "particle_position"]
   print ad["dinosaurs", "particle_position"]

Each of these three fields may have different sizes.  In order to enable
falling back on asking only for a field by the name, yt will use the most
recently requested field type for subsequent queries.  (By default, if no field
has been queried, it will look for the special field ``all``, which
concatenates all particle types.)  For example, if I were to then query for the
velocity:

.. code-block:: python

   print ad["particle_velocity"]

it would select ``dinosaurs`` as the field type.

The same operations work for fluid and mesh fields.  As an example, in some
cosmology simulations, we may want to examine the mass of particles in a region
versus the mass of gas.  We can do so by examining the special "deposit" field
types (described below) versus the gas fields:

.. code-block:: python

   print ad["deposit", "dark_matter_density"] / ad["gas", "density"]

The ``deposit`` field type is a mesh field, so it will have the same shape as
the gas density.  If we weren't using ``deposit``, and instead directly
querying a particle field, this *wouldn't* work, as they are different shapes.
This is the primary difference, in practice, between mesh and particle fields
-- they will be different shapes and so cannot be directly compared without
translating one to the other, typically through a "deposition" or "smoothing"
step.

How are fields implemented?
+++++++++++++++++++++++++++

There are two classes of fields in yt.  The first are those fields that exist
external to yt, which are immutable and can be queried -- most commonly, these
are fields that exist on disk.  These will often be returned in units that are
not in a known, external unit system (except possibly by design, on the part of
the code that wrote the data), and yt will take every effort possible to use
the names by which they are referred to by the data producer.  The default
field type for mesh fields that are "on-disk" is the name of the code frontend.
(For example, ``art``, ``enzo``, ``pyne``, and so on.) The default name for
particle fields, if they do not have a particle type affiliated with them, is
``io``.

The second class of field is the "derived field."  These are fields that are
functionally defined, either *ab initio* or as a transformation or combination
of other fields.  For example, when dealing with simulation codes, often the
fields that are evolved and output to disk are not the fields that are the most
relevant to researchers.  Rather than examining the internal gas energy, it is
more convenient to think of the temperature.  By applying one or multiple
functions to on-disk quantities, yt can construct new derived fields from them.
Derived fields do not always have to relate to the data found on disk; special
fields such as ``x``, ``y``, ``phi`` and ``dz`` all relate exclusively to the
geometry of the mesh, and provide information about the mesh that can be used
elsewhere for further transformations.

For more information, see :ref:`creating-derived-fields`.

There is a third, borderline class of field in yt, as well.  This is the
"alias" type, where a field on disk (for example, ``Density``) is aliased into
an internal yt-name (for example, ``density``).  The aliasing process allows
universally-defined derived fields to take advantage of internal names, and it
also provides an easy way to address what units something should be returned
in.  If an aliased field is requested (and aliased fields will always be
lowercase, with underscores separating words) it will be returned in CGS units
(future versions will enable global defaults to be set for MKS and other unit
systems), whereas if the underlying field is requested, it will not undergo any
unit conversions from its natural units.  (This rule is occasionally violated
for fields which are mesh-dependent, specifically particle masses in some
cosmology codes.)

Field types known to yt
+++++++++++++++++++++++

yt knows of a few different field types, by default.

 * ``index`` - this field type refers to characteristics of the mesh, whether
   that mesh is defined by the simulation or internally by an octree indexing
   of particle data.  A few handy fields are ``x``, ``y``, ``z``, ``theta``,
   ``phi``, ``radius``, ``dx``, ``dy``, ``dz`` and so on.
 * ``gas`` - this is the usual default for simulation frontends for fluid
   types.
 * ``all`` - this is a special particle field type that represents a
   concatenation of all particle field types.
 * ``deposit`` - this field type refers to the deposition of particles
   (discrete data) onto a mesh, typically to compute smoothing kernels, local
   density estimates, counts, and the like.
 * ``io`` - if a data frontend does not have a set of particle types, this will
   be the default for particle types.
 * frontend-name - mesh or fluid fields that exist on-disk default to having
   the name of the frontend as their type name. (i.e., ``enzo``, ``flash``,
   ``pyne`` and so on.)
 * particle type - if the particle types in the file are affiliated with names
   (rather than just ``io``) they will be available as field types.
   Additionally, any particle unions or filters will be accessible as field
   types.

Field Plugins
+++++++++++++

Derived fields are organized via plugins.  Inside yt are a number of field
plugins, which take information about fields in a dataset and then construct
derived fields on top of them.  This allows them to take into account
variations in naming system, units, data representations, and most importantly,
allows only the fields that are relevant to be added.  This system will be
expanded in future versions to enable much deeper semantic awareness of the
data types being analyzed by yt.

The field plugin system works in this order:

 * Available, inherent fields are identified by yt
 * The list of enabled field plugins is iterated over.  Each is called, and new
   derived fields are added as relevant.
 * Any fields which are not available, or which throw errors, are discarded.
 * Remaining fields are added to the list of derived fields available for a
   dataset
 * Dependencies for every derived field are identified, to enable data
   preloading

Field plugins can be loaded dynamically, although at present this is not
particularly useful.  Plans for extending field plugins to dynamically load, to
enable simple definition of common types (gradient, divergence, etc), and to
more verbosely describe available fields, have been put in place for future
versions.

The field plugins currently available include:

 * Angular momentum fields for particles and fluids
 * Astrophysical fields, such as those related to cosmology
 * Vector fields for fluid fields, such as gradients and divergences
 * Particle vector fields
 * Magnetic field-related fields
 * Species fields, such as for chemistry species (yt can recognize the entire
   periodic table in field names and construct ionization fields as need be)

What fields are available?
++++++++++++++++++++++++++

.. include reference here once it's done

The full list of fields available for a dataset can be found as 
the attribute ``field_list`` for native, on-disk fields and ``derived_field_list``
for derived fields (``derived_field_list`` is a superset of ``field_list``).
You can view these lists by examining a dataset like this:

.. code-block:: python

   ds = yt.load("my_data")
   print ds.field_list
   print ds.derived_field_list

By using the ``field_info()`` class, one can access information about a given
field, like its default units or the source code for it.  

.. code-block:: python

   ds = yt.load("my_data")
   ds.index
   print ds.field_info["gas", "pressure"].get_units()
   print ds.field_info["gas", "pressure"].get_source()

Particle Fields
---------------

Naturally, particle fields contain properties of particles rather than
grid cells.  Many of these fields have corresponding grid fields that
can be populated by "depositing" the particle values onto a yt grid.

General Particle Fields
+++++++++++++++++++++++

Every particle will contain both a ``particle_position`` and ``particle_velocity``
that tracks the position and velocity (respectively) in code units.


SPH Fields
++++++++++

For gas particles from SPH simulations, each particle will typically carry
a field for the smoothing length ``h``, which is roughly equivalent to 
``(m/\rho)^{1/3}``, where ``m`` and ``rho`` are the particle mass and density 
respectively.  This can be useful for doing neighbour finding.

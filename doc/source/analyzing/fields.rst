.. _fields:

Fields in yt
============

Fields are spatially-dependent quantities associated with a parent dataset.
Examples of fields are gas density, gas temperature, particle mass, etc.
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
changed to take an optional *field type* in addition to the *field name* of
the form ('*field type*', '*field name*').

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
---------------------------

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
"alias" type, where a field on disk (for example, (frontend, ``Density``)) is 
aliased into an internal yt-name (for example, (``gas``, ``density``)).  The 
aliasing process allows universally-defined derived fields to take advantage of 
internal names, and it also provides an easy way to address what units something 
should be returned in.  If an aliased field is requested (and aliased fields 
will always be lowercase, with underscores separating words) it will be returned 
in CGS units (future versions will enable global defaults to be set for MKS and 
other unit systems), whereas if the frontend-specific field is requested, it 
will not undergo any unit conversions from its natural units.  (This rule is 
occasionally violated for fields which are mesh-dependent, specifically particle 
masses in some cosmology codes.)

.. _known-field-types:

Field types known to yt
-----------------------

Recall that fields are formally accessed in two parts: ('*field type*', 
'*field name*').  Here we describe the different field types you will encounter:

* frontend-name -- Mesh or fluid fields that exist on-disk default to having
  the name of the frontend as their type name (e.g., ``enzo``, ``flash``,
  ``pyne`` and so on).  The units of these types are whatever units are
  designated by the source frontend when it writes the data.
* ``index`` -- This field type refers to characteristics of the mesh, whether
  that mesh is defined by the simulation or internally by an octree indexing
  of particle data.  A few handy fields are ``x``, ``y``, ``z``, ``theta``,
  ``phi``, ``radius``, ``dx``, ``dy``, ``dz`` and so on.  Default units
  are in CGS.
* ``gas`` -- This is the usual default for simulation frontends for fluid
  types.  These fields are typically aliased to the frontend-specific mesh
  fields for grid-based codes or to the deposit fields for particle-based
  codes.  Default units are in CGS.
* particle type -- These are particle fields that exist on-disk as written 
  by individual frontends.  If the frontend designates names for these particles
  (i.e. particle type) those names are the field types. 
  Additionally, any particle unions or filters will be accessible as field
  types.  Examples of particle types are ``Stars``, ``DM``, ``io``, etc.  
  Like the front-end specific mesh or fluid fields, the units of these fields
  are whatever was designated by the source frontend when written to disk.
* ``io`` -- If a data frontend does not have a set of multiple particle types, 
  this is the default for all particles.
* ``all`` -- This is a special particle field type that represents a
  concatenation of all particle field types using :ref:`particle-unions`.
* ``deposit`` -- This field type refers to the deposition of particles
  (discrete data) onto a mesh, typically to compute smoothing kernels, local
  density estimates, counts, and the like.  See :ref:`deposited-particle-fields` 
  for more information.

While it is best to be explicit access fields by their full names 
(i.e. ('*field type*', '*field name*')), yt provides an abbreviated 
interface for accessing common fields (i.e. '*field name*').  In the abbreviated
case, yt will assume you want the last *field type* accessed.  If you
haven't previously accessed a *field type*, it will default to *field type* = 
``'all'`` in the case of particle fields and *field type* = ``'gas'`` in the 
case of mesh fields.

Field Plugins
-------------

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
--------------------------

We provide a full list of fields that yt recognizes by default at 
:ref:`field-list`.  If you want to create additional custom derived fields, 
see :ref:`creating-derived-fields`.

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
grid cells.  By examining the particle field in detail, you can see that 
each element of the field array represents a single particle, whereas in mesh 
fields each element represents a single mesh cell.  This means that for the
most part, operations cannot operate on both particle fields and mesh fields
simultaneously in the same way, like filters (see :ref:`filtering-data`).
However, many of the particle fields have corresponding mesh fields that
can be populated by "depositing" the particle values onto a yt grid as 
described below.

.. _field_parameters:

Field Parameters
----------------

Certain fields require external information in order to be calculated.  For 
example, the radius field has to be defined based on some point of reference 
and the radial velocity field needs to know the bulk velocity of the data object 
so that it can be subtracted.  This information is passed into a field function 
by setting field parameters, which are user-specified data that can be associated 
with a data object.  The 
:meth:`~yt.data_objects.data_containers.YTDataContainer.set_field_parameter` 
and 
:meth:`~yt.data_objects.data_containers.YTDataContainer.get_field_parameter` 
functions are 
used to set and retrieve field parameter values for a given data object.  In the 
cases above, the field parameters are ``center`` and ``bulk_velocity`` respectively -- 
the two most commonly used field parameters.

.. code-block:: python

   ds = yt.load("my_data")
   ad = ds.all_data()

   ad.set_field_parameter("wickets", 13)

   print ad.get_field_parameter("wickets")

If a field parameter is not set, ``get_field_parameter`` will return None.  
Within a field function, these can then be retrieved and used in the same way.

.. code-block:: python

   def _wicket_density(field, data):
       n_wickets = data.get_field_parameter("wickets")
       if n_wickets is None:
           # use a default if unset
           n_wickets = 88
       return data["gas", "density"] * n_wickets

For a practical application of this, see :ref:`cookbook-radial-velocity`.

General Particle Fields
-----------------------

Every particle will contain both a ``particle_position`` and ``particle_velocity``
that tracks the position and velocity (respectively) in code units.

.. _deposited-particle-fields:

Deposited Particle Fields
-------------------------

In order to turn particle (discrete) fields into fields that are deposited in
some regular, space-filling way (even if that space is empty, it is defined
everywhere) yt provides mechanisms for depositing particles onto a mesh.  These
are in the special field-type space ``deposit``, and are typically of the form
``("deposit", "particletype_depositiontype")`` where ``depositiontype`` is the
mechanism by which the field is deposited, and ``particletype`` is the particle
type of the particles being deposited.  If you are attempting to examine the
cloud-in-cell (``cic``) deposition of the ``all`` particle type, you would
access the field ``("deposit", "all_cic")``.

yt defines a few particular types of deposition internally, and creating new
ones can be done by modifying the files ``yt/geometry/particle_deposit.pyx``
and ``yt/fields/particle_fields.py``, although that is an advanced topic
somewhat outside the scope of this section.  The default deposition types
available are:

* ``count`` - this field counts the total number of particles of a given type
  in a given mesh zone.  Note that because, in general, the mesh for particle
  datasets is defined by the number of particles in a region, this may not be
  the most useful metric.  This may be made more useful by depositing particle
  data onto an :ref:`arbitrary-grid`.
* ``density`` - this field takes the total sum of ``particle_mass`` in a given
  mesh field and divides by the volume.
* ``mass`` - this field takes the total sum of ``particle_mass`` in each mesh
  zone.
* ``cic`` - this field performs cloud-in-cell interpolation (see `Section 2.2
  <http://ta.twi.tudelft.nl/dv/users/Lemmens/MThesis.TTH/chapter4.html>`_ for more
  information) of the density of particles in a given mesh zone.
* ``smoothed`` - this is a special deposition type.  See discussion below for
  more information, in :ref:`sph-fields`.

.. _sph-fields:

SPH Fields
----------

For gas particles from SPH simulations, each particle will typically carry
a field for the smoothing length ``h``, which is roughly equivalent to 
``(m/\rho)^{1/3}``, where ``m`` and ``rho`` are the particle mass and density 
respectively.  This can be useful for doing neighbour finding.

As a note, SPH fields are special cases of the "deposited" particle fields.
They contain an additional piece of information about what is being examined,
and any fields that are recognized as being identical to intrinsic yt fields
will be aliased.  For example, in a Gadget dataset, the smoothed density of
``Gas`` particles will be aliased to the mesh field ``("gas", "density")`` so
that operations conducted on the mesh field ``density`` (which are frequent
occurrences) will operate on the smoothed gas density from the SPH particles.

The special deposition types based on smoothing (``smoothed``) are defined in
the file ``yt/geometry/particle_smooth.pyx``, and they require non-local
operations defined on a variable number of neighbors.  The default smoothing
type utilizes a cubic spline kernel and uses 64 nearest neighbors, providing a
volume-normalized smoothing.  Other types are possible, and yt provides
functionality for many different types of non-local correlation between
particles.  (For instance, a friends-of-friends grouper has been built on this
same infrastructure.)

Every particle field on a smoothed particle type is the source for a smoothed
field; this is not always useful, but it errs on the side of extra fields,
rather than too few fields.  (For instance, it may be unlikely that the
smoothed angular momentum field will be useful.)  The naming scheme is an
extension of the scheme described in :ref:`deposited-particle-fields`, and is
defined as such: ``("deposit", "particletype_smoothed_fieldname")``, where 
``fieldname`` is the name of the field being smoothed.  For example, smoothed
``Temperature`` of the ``Gas`` particle type would be ``("deposit",
"Gas_smoothed_Temperature")``, which in most cases would be aliased to the
field ``("gas", "temperature")`` for convenience.

Computing the Nth Nearest Neighbor
----------------------------------

One particularly useful field that can be created is that of the distance to
the Nth-nearest neighbor.  This field can then be used as input to smoothing
operations, in the case when a particular particle type does not have an
associated smoothing length or other length estimate.

yt defines this field as a plugin, and it can be added like so:

.. code-block:: python

   import yt
   from yt.fields.particle_fields import \
     add_nearest_neighbor_field

   ds = yt.load("snapshot_033/snap_033.0.hdf5")
   fn, = add_nearest_neighbor_field("all", "particle_position", ds)

   dd = ds.all_data()
   print dd[fn]

Note that ``fn`` here is the "field name" that yt adds.  It will be of the form
``(ptype, nearest_neighbor_distance_NN)`` where ``NN`` is the integer.  By
default this is 64, but it can be supplied as the final argument to
``add_nearest_neighbor_field``.  For the example above, it would be
``nearest_neighbor_64``.

This can then be used as input to the function
``add_volume_weighted_smoothed_field``, which can enable smoothing particle
types that would normally not be smoothed.

Gridded Data Format
===================

This is a pre-release of version 1.0 of this format.  Lots of formats have come
before, but this one is simple and will work with yt; the idea is to create an
import and export function in yt that will read this, so that other codes (such
as ZEUS-MP) can export directly to it or convert their data to it, and so that
yt can export to it from any format it recognizes and reads.

Caveats and Notes
-----------------

#. We avoid having many attributes on many nodes, as access can be quite slow
#. Cartesian data only for now
#. All grids must have the same number of ghost zones.
#. If “/grid_parent” does not exist, parentage relationships will be
   reconstructed and assumed to allow multiple grids
#. No parentage can skip levels
#. All grids are at the same time
#. This format is designed for single-fluid calculations (with color fields)
   but it should be viewed as extensible to multiple-fluids.
#. All fluid quantities are assumed to be in every grid, filling every zone.  Inside
   a given grid, for a given particle type, all the affiliated fields must be the
   same length.  (i.e., dark matter's velocity must be the same in all dimensions.)
#. Everything is in a single file; for extremely large datasets, the user may
   utilize HDF5 external links to link to files other than the primary.  (This
   enables, for instance, Enzo datasets to have only a thin wrapper that creates
   this format.)
#. All fluid fields in this version of the format are assumed to have the
   dimensionality of the grid they reside in plus any ghost zones, plus any
   additionaly dimensionality required by the staggering property.
#. Particles may have dataspaces affiliated with them.  (See Enzo's
   OutputParticleTypeGrouping for more information.)  This enables a light
   wrapper around data formats with interspersed particle types.
#. Boundary conditions are very simply specified -- future revisions
   will feature more complicated and rich specifications for the boundary.

Furthermore, we make a distinction between fluid quantities and particle
quantities.  Particles remain affiliated with grid nodes.  Positions of
particles are global, but this will change with future versions of this
document.

Format Declaration
------------------

The file type is HDF5.  We require version 1.8 or greater.  At the root level,
this group must exist: ::

   /gridded_data_format

This must contain the (float) attribute ``format_version``.  This document
describes version 1.0.  Optional attributes may exist:

``data_software``
   string, references the application creating the file, not the
   author of the data
``data_software_version``
   string, should reference a unique version number
``data_author``
   string, references the person or persons who created the data,
   should include an email address
``data_comment``
   string, anything about the data

Top Level Nodes
---------------

At least five top-level groups must exist, although some may be empty. ::

   /gridded_data_format
   /data
   /simulation_parameters
   /field_types
   /particle_types

Additionally, the grid structure elements must exist.  The 0-indexed index into this array
defines a unique "Grid ID".

``/grid_left_index``
   (int64, Nx3): global, relative to current level, and only the active region
``/grid_dimensions``
   (int64, Nx3): only the active regions
``/grid_level``
   (int64, N): level, indexed by zero
``/grid_particle_count``
   (int64, N): total number of particles.  (May change in subsequent versions.)
``/grid_parent_id``
   (int64, N): optional, may only reference a single parent

Grid Fields
-----------

Underneath ``/data/`` there must be entries for every grid, of the format
``/data/grid_%010i``.  These grids need no attributes, and underneath them
datasets live.

Fluid Fields
++++++++++++

For every grid we then define ``/data/grid_%010i/%(field)s``.

Where ``%(field)s`` draws from all of the fields defined.  We define no
standard for which fields must be present, only the names and units.  Units
should always be ''proper'' cgs (or conversion factors should be supplied, below), and
field names should be drawn from this list, with these names.  Not all fields
must be represented.  Field must extend beyond the active region if ghost zones
are included.  All pre-defined fields are assumed to be cell-centered unless this
is overridden in ``field_types``.

  * ``density`` (g/cc)
  * ``temperature`` (K)
  * ``specific_thermal_energy`` (erg/g)
  * ``specific_energy`` (erg/g, includes kinetic and magnetic)
  * ``magnetic_energy`` (erg/g)
  * ``velocity_x`` (cm/s)
  * ``velocity_y`` (cm/s)
  * ``velocity_z`` (cm/s)
  * ``species_density_%s`` (g/cc) where %s is the species name including ionization
    state, such as H2I, HI, HII, CO, "elec" for electron
  * ``mag_field_x``
  * ``mag_field_y``
  * ``mag_field_z``

Particle Fields
+++++++++++++++

Particles are more expensive to sort and identify based on "type" -- for
instance, dark matter versus star particles.  The particles should be separated
based on type, under the group ``/data/grid_%010i/particles/``.

The particles group will have sub-groups, each of which will be named after the
type of particle it represents.  We only specify "dark_matter" as a type;
anything else must be specified as described below.

Each node, for instance ``/data/grid_%010i/particles/dark_matter/``, must
contain the following fields:

  * ``mass`` (g)
  * ``id``
  * ``position_x`` (in physical units)
  * ``position_y`` (in physical units)
  * ``position_z`` (in physical units)
  * ``velocity_x`` (cm/s)
  * ``velocity_y`` (cm/s)
  * ``velocity_z`` (cm/s)
  * ``dataspace`` (optional) an HDF5 dataspace to be used when opening
    all affiliated fields.   If this is to be used, it must be appropriately set in
    the particle type definition.  This is of type ``H5T_STD_REF_DSETREG``.
    (See Enzo's OutputParticleTypeGrouping for an example.)

Additional Fields
+++++++++++++++++

Any additional fields from the data can be added, but must have a corresponding
entry in the root field table (described below.)  The naming scheme is to be as
explicit as possible, with units in cgs (or a conversion factor to the standard
cgs unit, in the field table.)

Attribute Table
---------------

In the root node, we define several groups which contain attributes.

Simulation Parameters
+++++++++++++++++++++

These attributes will all be associated with ``/simulation_parameters``.

``refine_by``
   relative global refinement
``dimensionality``
   1-, 2- or 3-D data
``domain_dimensions``
   dimensions in the top grid
``current_time``
   current time in simulation, in seconds, from “start” of simulation
``domain_left_edge``
   the left edge of the domain, in cm
``domain_right_edge``
   the right edge of the domain, in cm
``unique_identifier``
   regarded as a string, but can be anything
``cosmological_simulation``
   0 or 1
``num_ghost_zones``
   integer
``field_ordering``
   integer: 0 for C, 1 for Fortran
``boundary_conditions``
   integer (6): 0 for periodic, 1 for mirrored, 2 for outflow.  Needs one for each face
   of the cube.  Any past the dimensionality should be set to -1.  The order of specification
   goes left in 0th dimension, right in 0th dimension, left in 1st dimension, right in 1st dimensions,
   left in 2nd dimension, right in 2nd dimension.  Note also that yt does not currently support non-periodic
   boundary conditions, and that the assumption of periodicity shows up primarily in plots and
   covering grids.

Optionally, attributes for cosmological simulations can be provided, if
cosmological_simulation above is set to 1:

  * current_redshift
  * omega_matter (at z=0)
  * omega_lambda (at z=0)
  * hubble_constant (h100)

Fluid Field Attributes
++++++++++++++++++++++

Every field that is included that is not both in CGS already and in the list
above requires parameters.  If a field is in the above list but is not in CGS,
only the field_to_cgs attribute is necessary.  These will be stored under
``/field_types`` and each must possess the following attributes:

``field_name``
   a string that will be used to describe the field; can contain spaces.
``field_to_cgs``
   a float that will be used to convert the field to cgs units, if necessary.
   Set to 1.0 if no conversion necessary.  Note that if non-CGS units are desired
   this field should simply be viewed as the value by which field values are
   multiplied to get to some internally consistent unit system.
``field_units``
   a string that names the units.
``staggering``
   an integer: 0 for cell-centered, 1 for face-centered, 2 for vertex-centered.
   Non-cellcentered data will be linearly-interpolated; more complicated
   reconstruction will be defined in a future version of this standard; for 1.0
   we only allow for simple definitions.

Particle Types
++++++++++++++

Every particle type that is not recognized (i.e., all non-Dark Matter types)
needs to have an entry under ``/particle_types``.  Each entry must possess the
following attributes:

``particle_type_name``
   a string that will be used to describe the field; can contain spaces.
``particle_use_dataspace``
   (optional) if 1, the dataspace (see particle field definition above) will be used
   for all particle fields for this type of particle.  Useful if a given type of particle
   is embedded inside a larger list of different types of particle.
``particle_type_num``
   an integer giving the total number of particles of this type.

For instance, to define a particle of type ``accreting_black_hole``, the file
must contain ``/particle_types/accreting_black_hole``, with the
``particle_type_name`` attribute of "Accreting Black Hole".

Particle Field Attributes
+++++++++++++++++++++++++

Every particle type that contains a new field (for instance, ``accretion_rate``)
needs to have an entry under ``/particle_types/{particle_type_name}/{field_name}``
containing the following attributes:

``field_name``
   a string that will be used to describe the field; can contain spaces.
``field_to_cgs``
   a float that will be used to convert the field to cgs units, if necessary.
   Set to 1.0 if no conversion necessary.
``field_units``
   a string that names the units.

Role of YT
----------

yt will provide a reader for this data, so that any data in this format can be
used by the code.  Additionally, the names and specifications in this code
reflect the internal yt data structures.

yt will also provide a writer for this data, which will operate on any existing
data format.  Provided that a simulation code can read this data, this will
enable cross-platform comparison.  Furthermore, any external piece of software
(i.e., Stranger) that implements reading this format will be able to read any
format of data tha yt understands.

Example File
------------

An example file constructed from the ``RD0005-mine`` dataset is available
at http://yt.enzotools.org/files/RD0005.gdf .  It is not yet a complete
conversion, but it is a working proof of concept.  Readers and writers are
forthcoming.

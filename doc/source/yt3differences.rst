.. _yt3differences:

What's New and Different in yt 3.0?
===================================

If you are new to yt, welcome!  If you're coming to yt 3.0 from an older
version, however, there may be a few things in this version that are different
than what you are used to.  We have tried to build compatibility layers to
minimize disruption to existing scripts, but necessarily things will be
different in some ways.

Additionally, this list is current as of the latest alpha release.  Several
more changes are planned that are considerably more disruptive: these include
unit handling, index/index handling, and field naming.

.. warning:: This document covers *current* API changes.  As API changes occur
             it will be updated.

Cheat Sheet
-----------

Here's a quick reference for how to update your code to work with yt-3.0.

  * We have reworked yt's import system so that most commonly-used yt functions
    and classes live in the top-level ``yt`` namespace. That means you can now
    import yt with ``import yt``, load a dataset with ``ds = yt.load``
    and create a plot with ``yt.SlicePlot``.  See :ref:`api-reference` for a full
    API listing.
  * Fields and metadata for data objects and datasets now have units.  The unit
    system keeps you from making weird things like ``ergs`` + ``g`` and can
    handle things like ``g`` + ``kg`` or ``kg*m/s**2 == Newton``.  See
    :ref:`units` for more information.
  * Previously, yt would use "Enzo-isms" for field names. We now very
    specifically define fields as lowercase with underscores.  For instance,
    what used to be ``VelocityMagnitude`` would now be ``velocity_magnitude``.
    Axis names are now at the *end* of field names, not the beginning.
    ``x-velocity`` is now ``velocity_x``.
  * Fields can be accessed by a name, but are named internally as ``(fluid_type,
    fluid_name)``.
  * Mesh fields on-disk will be in code units, and will be named ``(code_name,
    FieldName)``.
  * Particle fields on-disk will also be in code units, and will be named
    ``(particle_type, FieldName)``.  If there is only one particle type in the
    output file, the particle type for all particles will be ``io``.
  * Previously, yt would capture command line arguments when being imported.
    This no longer happens.  As a side effect, it is no longer necessary to
    specify ``--parallel`` at the command line when running a parallel 
    computation. Use ``yt.enable_parallelism()`` instead.  See 
    :ref:`parallel-computation` for more detail.
  * Any derived quantities that *always* returned lists (like ``Extrema``,
    which would return a list even if you only ask for one field) now only
    returns a single result if you only ask for one field.  Results for particle
    and mesh fields will be returned separately.
  * Derived quantities can now be accessed via a function that hangs off of the
    ``quantities`` atribute of data objects. Instead of
    ``dd.quantities['TotalMass']``, you can now use
    ``dd.quantities.total_mass()`` to do the same thing. All derived quantities
    can be accessed via a function that hangs off of the `quantities` attribute
    of data objects.

Cool New Things
---------------

Lots of new things have been added in yt 3.0!  Below we summarize a handful of
these.

Octrees
+++++++

Octree datasets such as RAMSES, ART and ARTIO are now supported -- without any
regridding!  We have a native, lightweight octree indexing system.

Particle Codes and SPH
++++++++++++++++++++++

yt 3.0 features particle selection, smoothing, and deposition.  This utilizes a
combination of coarse-grained indexing and octree indexing for particles.

Irregular Grids
+++++++++++++++

MOAB Hex8 format is supported, and non-regular grids can be added relatively
easily.

Particle Deposition
+++++++++++++++++++

In yt-3.0, we provide mechanisms for describing and creating fields generated
by depositing particles into one or a handful of zones.  This could include
deposited mass or density, average values, and the like.  For instance, the
total stellar mass in some region can be deposited and averaged.

Particle Filters and Unions
+++++++++++++++++++++++++++

Throughout yt, the notion of "particle types" has been more deeply embedded.
These particle types can be dynamically defined at runtime, for instance by
taking a filter of a given type or the union of several different types.  This
might be, for instance, defining a new type called ``young_stars`` that is a
filtering of ``star_age`` to be fewer than a given threshold, or ``fast`` that
filters based on the velocity of a particle.  Unions could be the joining of
multiple types of particles -- the default union of which is ``all``,
representing all particle types in the simulation.

Units
+++++

yt now has a unit system.  This is one of the bigger features, and in essence it means
that you can convert units between anything.  In practice, it makes it much
easier to define fields and convert data between different unit systems. See
:ref:`units` for more information.

Non-Cartesian Coordinates
+++++++++++++++++++++++++

Preliminary support for non-cartesian coordinates has been added.  We expect
this to be considerably solidified and expanded in yt 3.1.

API Changes
-----------

These are the items that have already changed in *user-facing* API:

Field Naming
++++++++++++

.. warning:: Field naming is probably the single biggest change you will
             encounter in yt 3.0.

Fields can be accessed by their short names, but yt now has an explicit
mechanism of distinguishing between field types and particle types.  This is
expressed through a two-key description.  For example::

   my_object["gas", "density"]

will return the gas field density.  In this example "gas" is the field type and
"density" is the field name.  Field types are a bit like a namespace.  This
system extends to particle types as well.  By default you do *not* need to use
the field "type" key, but in case of ambiguity it will utilize the default value
in its place.  This should therefore be identical to::

   my_object["density"]

Units of Fields
+++++++++++++++

Fields now are all subclasses of NumPy arrays, the ``YTArray``, which carries
along with it units.  This means that if you want to manipulate fields, you
have to modify them in a unitful way.

Parameter Files are Now Datasets
++++++++++++++++++++++++++++++++

Wherever possible, we have attempted to replace the term "parameter file"
(i.e., ``pf``) with the term "dataset."  Future revisions will change most of
the ``pf`` atrributes of objects into ``ds`` or ``dataset`` attributes.

Hierarchy is Now Index
++++++++++++++++++++++

The hierarchy object (``pf.h``) is now referred to as an index (``ds.index``).
It is no longer necessary to directly refer to the ``index`` as often, since
data objects are now attached to the to the ``dataset`` object.  Before, you
would say ``ph.f.sphere()``, now you can say ``ds.sphere()``.

Field Info
++++++++++

In previous versions of yt, the ``dataset`` object (what we used to call a
parameter file) had a ``field_info`` attribute which was a dictionary leading to
derived field definitions.  At the present time, because of the field naming
changes (i.e., access-by-tuple) it is better to utilize the function
``_get_field_info`` than to directly access the ``field_info`` dictionary.  For
example::

   finfo = ds._get_field_info("gas", "density")

This function respects the special "field type" ``unknown`` and will search all
field types for the field name.

Projection Argument Order
+++++++++++++++++++++++++

Previously, projections were inconsistent with the other data objects.
(The API for Plot Windows is the same.)  The argument order is now ``field``
then ``axis``.

Field Parameters
++++++++++++++++

All data objects now accept an explicit list of ``field_parameters`` rather
than accepting ``kwargs`` and supplying them to field parameters.

Object Renaming
+++++++++++++++

Nearly all internal objects have been renamed.  Typically this means either
removing ``AMR`` from the prefix or replacing it with ``YT``.  All names of
objects remain the same for the purposes of selecting data and creating them;
i.e., ``sphere`` objects are still called ``sphere`` - you can access create one
via ``ds.sphere``.

Boolean Regions
+++++++++++++++

Boolean regions are not yet implemented in yt 3.0.

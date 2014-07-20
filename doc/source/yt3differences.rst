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

  * Importing yt is now as simple as ``import yt``.  The docs have been
    extensively updated to reflect this new style.  ``from yt.mods import *``
    still works, but we are discouraging its use going forward.
  * Fields can be accessed by a name, but are named internally as ``(fluid_type,
    fluid_name)``.
  * Fields on-disk will be in code units, and will be named ``(code_name,
    FieldName)``.
  * Previously, yt would use "Enzo-isms" for field names.  We now very
    specifically define fields as lowercase with underscores.  For instance,
    what used to be ``VelocityMagnitude`` would now be ``velocity_magnitude``.
  * Particles are either named by their type or default to the type ``io``.
  * Axis names are now at the *end* of field names, not the beginning.
    ``x-velocity`` is now ``velocity_x``.
  * Any derived quantities that *always* returned lists (like ``Extrema``,
    which would return a list even if you only ask for one field) now only
    return a single tuple if you only ask for one field.
  * Units can be tricky, and they try to keep you from making weird things like
    ``ergs`` + ``g``.  See :ref:`units` for more information.
  * Previously, yt would capture command line arguments when being imported.
    This no longer happens.  As a side effect, it is no longer necessary to
    specify ``--parallel`` at the command line when running a parallel 
    computation. Use ``yt.enable_parallelism()`` instead.  See 
    :ref:`parallel-computation` for more detail.

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

yt now has units.  This is one of the bigger features, and in essence it means
that you can convert units between anything.  See :ref:`units` for more
information.

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

will return the gas field density.  This extends to particle types as well.  By
default you do *not* need to use the field "type" key, but in case of ambiguity
it will utilize the default value in its place.  This should therefore be
identical to::

   my_object["density"]

Units of Fields
+++++++++++++++

Fields now are all subclasses of NumPy arrays, the ``YTArray``, which carries
along with it units.  This means that if you want to manipulate fields, you
have to modify them in a unitful way.

Field Info
++++++++++

In the past, the object ``ds`` (or ``ds``) had a ``field_info`` object which
was a dictionary leading to derived field definitions.  At the present time,
because of the field naming changes (i.e., access-by-tuple) it is better to
utilize the function ``_get_field_info`` than to directly access the
``field_info`` dictionary.  For example::

   finfo = ds._get_field_info("gas", "density")

This function respects the special "field type" ``unknown`` and will search all
field types for the field name.

Parameter Files are Now Datasets
++++++++++++++++++++++++++++++++

Wherever possible, we have attempted to replace the term "parameter file"
(i.e., ``ds``) with the term "dataset."  Future revisions will change most of
the ``ds`` atrributes of objects into ``ds`` or ``dataset`` attributes.

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
i.e., you will not need to change ``ds.sphere`` to something else.

Boolean Regions
+++++++++++++++

Boolean regions are not yet implemented in yt 3.0.

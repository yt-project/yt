.. _yt3differences:

What's New and Different in yt 3.0?
===================================

If you are new to yt, welcome!  If you're coming to yt 3.0 from an older
version, however, there may be a few things in this version that are different
than what you are used to.  We have tried to build compatibility layers to
minimize disruption to existing scripts, but necessarily things will be
different in some ways.

.. contents::
   :depth: 2
   :local:
   :backlinks: none

Updating to yt 3.0 from Old Versions (and going back)
-----------------------------------------------------

First off, you need to update your version of yt to yt 3.0.  If you're
installing yt for the first time, please visit :ref:`getting-and-installing-yt`.
If you already have a version of yt installed, you should just need one
command:

.. code-block:: bash

    $ yt update

This will update yt to the most recent version and rebuild the source base.
If you installed using the installer script, it will assure you have all of the
latest dependencies as well.  This step may take a few minutes.  To test
to make sure yt is running, try:

.. code-block:: bash

    $ yt --help

If you receive no errors, then you are ready to go.  If you have
an error, then consult :ref:`update-errors` for solutions.

If you want to switch back to an old version of yt (2.x), see
:ref:`switching-between-yt-versions`.

.. _transitioning-to-3.0:

Converting Old Scripts to Work with yt 3.0
------------------------------------------

After installing yt-3.0, you'll want to change your old scripts in a few key
ways.  After accounting for the changes described in the list below, try
running your script.  If it still fails, the callback failures in python are
fairly descriptive and it may be possible to deduce what remaining changes are
necessary.  If you continue to have trouble, please don't hesitate to
:ref:`request help <asking-for-help>`.

The list below is arranged in order of most important changes to least
important changes.

* **Replace** ``from yt.mods import *`` **with** ``import yt`` **and prepend yt
  classes and functions with** ``yt.``
  We have reworked yt's import system so that most commonly-used yt functions
  and classes live in the top-level yt namespace. That means you can now
  import yt with ``import yt``, load a dataset with ``ds = yt.load(filename)``
  and create a plot with ``yt.SlicePlot``.  See :ref:`api-reference` for a full
  API listing.  You can still import using ``from yt.mods import *`` to get a
  pylab-like experience.
* **Unit conversions are different**
  Fields and metadata for data objects and datasets now have units.  The unit
  system keeps you from making weird things like ``ergs`` + ``g`` and can
  handle things like ``g`` + ``kg`` or ``kg*m/s**2 == Newton``.  See
  :ref:`units` and :ref:`conversion-factors` for more information.
* **Change field names from CamelCase to lower_case_with_underscores**
  Previously, yt would use "Enzo-isms" for field names. We now very
  specifically define fields as lowercase with underscores.  For instance,
  what used to be ``VelocityMagnitude`` would now be ``velocity_magnitude``.
  Axis names are now at the *end* of field names, not the beginning.
  ``x-velocity`` is now ``velocity_x``.  For a full list of all of the fields,
  see :ref:`field-list`.
* **Full field names have two parts now**
  Fields can be accessed by a single name, but they are named internally as
  ``(field_type, field_name)`` for more explicit designation which can address
  particles, deposited fluid quantities, and more.  See :ref:`fields`.
* **Code-specific field names can be accessed by the name defined by the
  external code**
  Mesh fields that exist on-disk in an output file can be read in using whatever
  name is used by the output file.  On-disk fields are always returned in code
  units.  The full field name will be ``(code_name, field_name)``. See
  :ref:`field-list`.
* **Particle fields are now more obviously different than mesh fields**
  Particle fields on-disk will also be in code units, and will be named
  ``(particle_type, field_name)``.  If there is only one particle type in the
  output file, all particles will use ``io`` as the particle type. See
  :ref:`fields`.
* **Change** ``pf`` **to** ``ds``
  The objects we used to refer to as "parameter files" we now refer to as
  datasets.  Instead of ``pf``, we now suggest you use ``ds`` to refer to an
  object returned by ``yt.load``.
* **Remove any references to** ``pf.h`` **with** ``ds``
  You can now create data objects without referring to the hierarchy. Instead
  of ``pf.h.all_data()``, you can now say ``ds.all_data()``.  The hierarchy is
  still there, but it is now called the index: ``ds.index``.
* **Use** ``yt.enable_parallelism()`` **to make a script parallel-compatible**
  Command line arguments are only parsed when yt is imported using ``from
  yt.mods import *``. Since command line arguments are not parsed when using
  ``import yt``, it is no longer necessary to specify ``--parallel`` at the
  command line when running a parallel computation. Use
  ``yt.enable_parallelism()`` in your script instead.  See
  :ref:`parallel-computation` for more details.
* **Change your derived quantities to the new syntax**
  Derived quantities have been reworked.  You can now do
  ``dd.quantities.total_mass()`` instead of ``dd.quantities['TotalMass']()``.
  See :ref:`derived-quantities`.
* **Change your method of accessing the** ``grids`` **attribute**
  The ``grids`` attribute of data objects no longer exists.  To get this
  information, you have to use spatial chunking and then access them.  See
  :ref:`here <grid-chunking>` for an example.  For datasets that use grid
  hierarchies, you can also access the grids for the entire dataset via
  ``ds.index.grids``.  This attribute is not defined for particle or octree
  datasets.

Cool New Things
---------------

Lots of new things have been added in yt 3.0!  Below we summarize a handful of
these.

Lots of New Codes are Supported
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because of the additions of **Octrees**, **Particle Deposition**,
and **Irregular Grids**, we now support a bunch more codes.  See
:ref:`code-support` for more information.

Octrees
^^^^^^^

Octree datasets such as RAMSES, ART and ARTIO are now supported -- without any
regridding!  We have a native, lightweight octree indexing system.

Irregular Grids
^^^^^^^^^^^^^^^

MOAB Hex8 format is supported, and non-regular grids can be added relatively
easily.

Better Particle Support
^^^^^^^^^^^^^^^^^^^^^^^

Particle Codes and SPH
""""""""""""""""""""""

yt 3.0 features particle selection, smoothing, and deposition.  This utilizes a
combination of coarse-grained indexing and octree indexing for particles.

Particle Deposition
"""""""""""""""""""

In yt-3.0, we provide mechanisms for describing and creating fields generated
by depositing particles into one or a handful of zones.  This could include
deposited mass or density, average values, and the like.  For instance, the
total stellar mass in some region can be deposited and averaged.

Particle Filters and Unions
"""""""""""""""""""""""""""

Throughout yt, the notion of "particle types" has been more deeply embedded.
These particle types can be dynamically defined at runtime, for instance by
taking a filter of a given type or the union of several different types.  This
might be, for instance, defining a new type called ``young_stars`` that is a
filtering of ``star_age`` to be fewer than a given threshold, or ``fast`` that
filters based on the velocity of a particle.  Unions could be the joining of
multiple types of particles -- the default union of which is ``all``,
representing all particle types in the simulation.

Units
^^^^^

yt now has a unit system.  This is one of the bigger features, and in essence it means
that you can convert units between anything.  In practice, it makes it much
easier to define fields and convert data between different unit systems. See
:ref:`units` for more information.

Non-Cartesian Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^

Preliminary support for non-cartesian coordinates has been added.  We expect
this to be considerably solidified and expanded in yt 3.1.

Reworked Import System
^^^^^^^^^^^^^^^^^^^^^^

It's now possible to import all yt functionality using ``import yt``. Rather
than using ``from yt.mods import *``, we suggest using ``import yt`` in new
scripts.  Most commonly used yt functionality is attached to the ``yt`` module.
Load a dataset with ``yt.load()``, create a phase plot using ``yt.PhasePlot``,
and much more, see :ref:`the api docs <api-reference>` to learn more about what's
in the ``yt`` namespace, or just use tab completion in IPython: ``yt.<tab>``.

It's still possible to use ``from yt.mods import *`` to create an interactive
pylab-like experience.  Importing yt this way has several side effects, most
notably the command line arguments parsing and other startup tasks will run.

API Changes
-----------

These are the items that have already changed in *user-facing* API:

Field Naming
^^^^^^^^^^^^

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

To enable a compatibility layer, on the dataset you simply need to call the
method ``setup_deprecated_fields`` like so:

.. code-block:: python

   ds = yt.load("MyData")
   ds.setup_deprecated_fields()

This sets up aliases from the old names to the new.  See :ref:`fields` and
:ref:`field-list` for more information.

Units of Fields
^^^^^^^^^^^^^^^

Fields now are all subclasses of NumPy arrays, the ``YTArray``, which carries
along with it units.  This means that if you want to manipulate fields, you
have to modify them in a unitful way.  See :ref:`units`.

Parameter Files are Now Datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wherever possible, we have attempted to replace the term "parameter file"
(i.e., ``pf``) with the term "dataset."  In yt-3.0, all of
the ``pf`` attributes of objects are now ``ds`` or ``dataset`` attributes.

Hierarchy is Now Index
^^^^^^^^^^^^^^^^^^^^^^

The hierarchy object (``pf.h``) is now referred to as an index (``ds.index``).
It is no longer necessary to directly refer to the ``index`` as often, since
data objects are now attached to the to the ``dataset`` object.  Before, you
would say ``pf.h.sphere()``, now you can say ``ds.sphere()``.

New derived quantities interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Derived quantities can now be accessed via a function that hangs off of the
``quantities`` attribute of data objects. Instead of
``dd.quantities['TotalMass']()``, you can now use ``dd.quantities.total_mass()``
to do the same thing. All derived quantities can be accessed via a function that
hangs off of the ``quantities`` attribute of data objects.

Any derived quantities that *always* returned lists (like ``Extrema``, which
would return a list even if you only ask for one field) now only returns a
single result if you only ask for one field.  Results for particle and mesh
fields will also be returned separately.  See :ref:`derived-quantities` for more
information.


Field Info
^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^

Previously, projections were inconsistent with the other data objects.
(The API for Plot Windows is the same.)  The argument order is now ``field``
then ``axis`` as seen here:
:class:`~yt.data_objects.construction_data_containers.YTQuadTreeProj`.

Field Parameters
^^^^^^^^^^^^^^^^

All data objects now accept an explicit list of ``field_parameters`` rather
than accepting ``kwargs`` and supplying them to field parameters.  See
:ref:`field_parameters`.

Object Renaming
^^^^^^^^^^^^^^^

Nearly all internal objects have been renamed.  Typically this means either
removing ``AMR`` from the prefix or replacing it with ``YT``.  All names of
objects remain the same for the purposes of selecting data and creating them;
i.e., ``sphere`` objects are still called ``sphere`` - you can access or create one
via ``ds.sphere``.  For a detailed description and index see
:ref:`available-objects`.

Boolean Regions
^^^^^^^^^^^^^^^

Boolean regions are not yet implemented in yt 3.0.

.. _grid-chunking:

Grids
^^^^^

It used to be that one could get access to the grids that belonged to a data
object.  Because we no longer have just grid-based data in yt, this attribute
does not make sense.  If you need to determine which grids contribute to a
given object, you can either query the ``grid_indices`` field, or mandate
spatial chunking like so:

.. code-block:: python

   for chunk in obj.chunks([], "spatial"):
       for grid in chunk._current_chunk.objs:
           print(grid)

This will "spatially" chunk the ``obj`` object and print out all the grids
included.

Halo Catalogs
^^^^^^^^^^^^^

The ``Halo Profiler`` infrastructure has been fundamentally rewritten and now
exists using the ``Halo Catalog`` framework.  See :ref:`halo-analysis`.

Analysis Modules
^^^^^^^^^^^^^^^^

While we're trying to port over all of the old analysis modules, we have not
gotten all of them working in 3.0 yet.  The docs pages for those modules
not-yet-functioning are clearly marked.

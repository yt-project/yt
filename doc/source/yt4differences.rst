.. _yt4differences:

What's New and Different in yt 4.0?
===================================

If you are new to yt, welcome!  If you're coming to yt 4.0 from an older
version, however, there may be a few things in this version that are different
than what you are used to.  We have tried to build compatibility layers to
minimize disruption to existing scripts, but necessarily things will be
different in some ways.

.. contents::
   :depth: 2
   :local:
   :backlinks: none

Updating to yt 4.0 from Old Versions (and going back)
-----------------------------------------------------


.. _transitioning-to-4.0:

Converting Old Scripts to Work with yt 4.0
------------------------------------------

After installing yt-4.0, you’ll want to change your old scripts in a few key
ways. After accounting for the changes described in the list below, try
running your script. If it still fails, the Python tracebacks
should be fairly descriptive and it may be possible to deduce what remaining
changes are necessary. If you continue to have trouble, please don’t hesitate
to :ref:`request help <asking-for-help>`.

The list below is arranged in order of most to least important changes.

* **Fields should be specified as tuples not as strings**
  In the past, you could specify fields as strings like ``"density"``, but
  with the growth of yt and its many derived fields, there can be sometimes
  be overlapping field names (e.g., ``("gas", "density")`` and
  ``("PartType0", "density")``, where yt doesn't know which to use.  To remove
  any ambiguity, it is now strongly recommended to explicitly specify the full
  tuple form of all fields. Just search for all field accesses in your scripts,
  and replace strings with tuples (e.g. replace ``"a"``  with
  ``("gas", "a" )``).  There is a compatibility rule in yt-4.0 to allow strings
  to continue to work until yt-4.1, but you may get unexpected behavior.  Any
  field specifications that are ambiguous will throw an error in future
  versions of yt.  See our :ref:`fields`, and :ref:`available field list
  <available-fields>` documentation for more information.
* **Use Newer Versions of Python**
  The yt-4.0 release will be the final release of yt to support Python 3.6.
  Starting with yt-4.1, python 3.6 will no longer be supported, so please
  start using 3.7+ as soon as possible.
* **Particle-based datasets no longer accept n_ref and over_refine_factor**
  One of the major upgrades in yt-4 is native treatment of particle-based
  datasets.  This is in contrast to previous yt behavior which loaded particle-based
  datasets as octrees, which could then be treated like grid-based datasets.
  In order to define the octrees, users were required to specify ``n_ref``
  and ``over_refine_factor`` values at load time.  Please remove
  any reference to ``n_ref`` and ``over_refine_factor`` in your scripts.
* **Neutral ion fields changing format**
  In previous versions, neutral ion fields were specified as
  ``ELEMENT_number_density`` (e.g., ``H_number_density`` to represent H I
  number density).  This led to a lot of confusion, because some people assumed
  these fields were the total hydrogen density, not neutral hydrogen density.
  In yt-4.0, we have resolved this issue by explicitly calling total hydrogen
  number density ``H_nuclei_density`` and neutral hydrogen density
  ``H_p0_number_density`` (where ``p0`` refers to plus 0 charge).  This syntax
  follows the rule for other ions: H II = ``H_p1`` = ionized hydrogen.  Change
  your scripts accordingly.  See :ref:`species-fields` for more information.
* **Change in energy and momentum field names**
  Fields representing energy and momentum quantities are now given names which
  reflect their dimensionality. For example, the ``("gas", "kinetic_energy")``
  field was actually a field for kinetic energy density, and so it has been
  renamed to ``"gas", "kinetic_energy_density"``. The old name still exists
  as an alias as of yt v4.0.0, but it will be removed in yt v4.1.0. See
  :ref:`deprecated_field_names` below for more information.
  Other examples include ``"gas", "specific_thermal_energy"`` for thermal
  energy per unit mass, and ``("gas", "momentum_density_x")`` for the x-axis
  component of momentum density. See :ref:`efields` for more information.
* **Deprecated field names**
  Certain field names are deprecated within yt v4.0 and will be removed in
  yt v4.1. For example, ``("gas", "kinetic_energy")`` has been renamed to
  ``("gas", "kinetic_energy_density")``, though the former name has been added
  as an alias. Other fields, such as
  ``("gas", "cylindrical_tangential_velocity_absolute")``, are being removed
  entirely. When the deprecated field names are used for the first time in a
  session, a warning will be logged, so it is advisable to set
  your logging level to ``WARNING`` (``yt.set_log_level("error")``) at a
  minimum to catch these.  See :ref:`faq-log-level` for more information on
  setting your log level and :ref:`available-fields` to see all available
  fields.
* ``cmocean`` **colormaps need prefixing**
  yt used to automatically load and register external colormaps from the
  ``cmocean`` package unprefixed (e.g., ``set_cmap(FIELD, "balance")``.  This
  became unsustainable with the 3.4 release of Matplotlib, in which colormaps
  with colliding names raise errors. The fix is to explicitly import the
  ``cmocean`` module and prefix ``cmocean`` colormaps (like ``balance``) with
  ``cmo.`` (e.g., ``cmo.balance``).  Note that this solution works with any
  yt-supported version of Matplotlib, but is not backward compatible with
  earlier versions of yt.
* Position and velocity fields now default to using linear scaling in profiles
  and phase plots, whereas previously behavior was determined by whether the
  dataset was particle- or grid-based.  Efforts have been made to standardize
  the treatment of other fields in profile and phase plots for particle and
  grid datasets.

Important New Aliases
^^^^^^^^^^^^^^^^^^^^^

With the advent of supporting SPH data at the particle level instead of smoothing
onto an octree (see below), a new alias for both gas particle masses and cell masses
has been created: ``("gas", "mass")``, which aliases to ``("gas", "cell_mass")`` for
grid-based frontends and to the gas particle mass for SPH frontends. In a number of
places in yt, code that used ``("gas", "cell_mass")`` has been replaced by
``("gas", "mass")``. Since the latter is an alias for the former, old scripts which
use ``("gas", "cell_mass")`` should not break.

Cool New Things
---------------

Changes for Working with SPH Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In yt-3.0 most user-facing operations on SPH data are produced by interpolating
SPH data onto a volume-filling octree mesh. Historically this was easier to
implement When support for SPH data was added to yt as it allowed re-using a lot
of the existing infrastructure. This had some downsides because the octree was a
single, global object, the memory and CPU overhead of smoothing SPH data onto
the octree can be prohibitive on particle datasets produced by large
simulations. Constructing the octree during the initial indexing phase also
required each particle (albeit, in a 64-bit integer) to be present in memory
simultaneously for a sorting operation, which was memory prohibitive.
Visualizations of slices and projections produced by yt using the default
settings are somewhat blocky since by default we use a relatively coarse octree
to preserve memory.

In yt-4.0 this has all changed! Over the past two years, Nathan Goldbaum, Meagan
Lang and Matt Turk implemented a new approach for handling I/O of particle data,
based on storing compressed bitmaps containing Morton indices instead of an
in-memory octree. This new capability means that the global octree index is now
no longer necessary to enable I/O chunking and spatial indexing of particle data
in yt.

The new I/O method has opened up a new way of dealing with the particle data and
in particular, SPH data.

Scatter and Gather approach for SPH data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As mentioned, previously operations such as slice, projection and arbitrary
grids would smooth the particle data onto the global octree. As this is no
longer used, a different approach was required to visualize the SPH data. Using
SPLASH as inspiration, SPH smoothing pixelization operations were created using
smooting operations via "scatter" and "gather" approaches. We estimate the
contributions of a particle to a single pixel by considering the point at the
centre of the pixel and using the standard SPH smoothing formula. The heavy
lifting in these functions is undertaken by cython functions.

It is now possible to generate slice plots, projection plots, covering grids and
arbitrary grids of smoothed quanitities using these operations. The following
code demonstrates how this could be achieved. The following would use the scatter
method:

.. code-block:: python

    import yt

    ds = yt.load("snapshot_033/snap_033.0.hdf5")

    plot = yt.SlicePlot(ds, 2, ("gas", "density"))
    plot.save()

    plot = yt.ProjectionPlot(ds, 2, ("gas", "density"))
    plot.save()

    arbitrary_grid = ds.arbitrary_grid([0.0, 0.0, 0.0], [25, 25, 25], dims=[16, 16, 16])
    ag_density = arbitrary_grid[("gas", "density")]

    covering_grid = ds.covering_grid(4, 0, 16)
    cg_density = covering_grid[("gas", "density")]

In the above example the ``covering_grid`` and the ``arbitrary_grid`` will return
the same data. In fact, these containers are very similar but provide a
slighlty different API.

The above code can be modified to use the gather approach by changing a global
setting for the dataset. This can be achieved with
``ds.sph_smoothing_style = "gather"``, so far, the gather approach is not
supported for projections.

The default behaviour for SPH interpolation is that the values are normalized
inline with Eq. 9 in `SPLASH, Price (2009) <https://arxiv.org/pdf/0709.0832.pdf>`_.
This can be disabled with ``ds.use_sph_normalization = False``. This will
disable the normalization for all future interpolations.

The gather approach requires finding nearest neighbors using the KDTree. The
first call will generate a KDTree for the entire dataset which will be stored in
a sidecar file. This will be loaded whenever neccesary.

Off-Axis Projection for SPH Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current ``OffAxisProjectionPlot`` class will now support SPH projection plots.

The following is a code example:

.. code-block:: python

    import yt

    ds = yt.load("Data/GadgetDiskGalaxy/snapshot_200.hdf5")

    smoothing_field = ("gas", "density")

    _, center = ds.find_max(smoothing_field)

    sp = ds.sphere(center, (10, "kpc"))

    normal_vector = sp.quantities.angular_momentum_vector()

    prj = yt.OffAxisProjectionPlot(ds, normal_vector, smoothing_field, center, (20, "kpc"))

    prj.save()

Smoothing Data onto an Octree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whilst the move away from the global octree is a promising one in terms of
perfomance and dealing with SPH data in a more intuitive manner, it does remove
a useful feature. We are aware that many users will have older scripts which take
advantage of the global octree.

As such, we have added support to smooth SPH data onto an octree when desired by
the users. The new octree is designed to give results consistent with those of
the previous octree, but the new octree takes advantage of the scatter and
gather machinery also added.

.. code-block:: python

    import numpy as np

    import yt

    ds = yt.load("GadgetDiskGalaxy/snapshot_200.hdf5")
    left = np.array([0, 0, 0], dtype="float64")
    right = np.array([64000, 64000, 64000], dtype="float64")

    # generate an octree
    octree = ds.octree(left, right, n_ref=64)

    # Scatter deposition is the default now, and thus this will print scatter
    print(octree.sph_smoothing_style)

    # the density will be calculated using SPH scatter
    density = octree[("PartType0", "density")]

    # this will return the x positions of the octs
    x = octree[("index", "x")]

The above code can be modified to use the gather approach by using
``ds.sph_smoothing_style = 'gather'`` before any field access. The octree just
uses the smoothing style and number of neighbors defined by the dataset.

The octree implementation is very simple. It uses a recursive algorithm to build
a ``depth-first`` which is consistent with the results from yt-3. Depth-first
search (DFS) means that tree starts refining at the root node (this is the
largest node which contains every particles) and refines as far as possible
along each branch before backtracking.

``yt.units`` Is Now a Wrapper for ``unyt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have extracted ``yt.units`` into ``unyt``, its own library that you can
install separately from yt from ``pypi`` and ``conda-forge``. You can find out
more about using ``unyt`` in `its documentation
<https://unyt.readthedocs.io/en/stable/>`_ and in `a paper in the Journal of
Open Source Software <http://joss.theoj.org/papers/10.21105/joss.00809>`_.

From the perspective of a user of yt, very little should change. While things in
``unyt`` have different names -- for example ``YTArray`` is now called
``unyt_array`` -- we have provided wrappers in ``yt.units`` so imports in your
old scripts should continue to work without issue. If you have any old scripts
that don't work due to issues with how yt is using ``unyt`` or units issues in
general please let us know by `filing an issue on GitHub
<https://github.com/yt-project/yt/issues/new>`_.

Moving ``unyt`` into its own library has made it much easier to add some cool
new features, which we detail below.

``ds.units``
~~~~~~~~~~~~

Each dataset now has a set of unit symbols and physical constants associated
with it, allowing easier customization and smoother interaction, especially in
workflows that need to use code units or cosmological units. The ``ds.units``
object has a large number of attributes corresponding to the names of units and
physical constants. All units known to the dataset will be available, including
custom units. In situations where you might have used ``ds.arr`` or ``ds.quan``
before, you can now safely use ``ds.units``:

   >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
   >>> u = ds.units
   >>> ad = ds.all_data()
   >>> data = ad['Enzo', 'Density']
   >>> data + 12*u.code_mass/u.code_length**3
   unyt_array([1.21784693e+01, 1.21789148e+01, 1.21788494e+01, ...,
               4.08936836e+04, 5.78006836e+04, 3.97766906e+05], 'code_mass/code_length**3')
   >>> data + .0001*u.mh/u.cm**3
   unyt_array([6.07964513e+01, 6.07968968e+01, 6.07968314e+01, ...,
               4.09423016e+04, 5.78493016e+04, 3.97815524e+05], 'code_mass/code_length**3')


Automatic Unit Simplification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often the results of an operation will result in a unit expression that can be
simplified by cancelling pairs of factors. Before yt 4.0, these pairs of factors
were only cancelled if the same unit appeared in both the numerator and
denominator of an expression. Now, all pairs of factors have have inverse
dimensions are cancelled, and the appropriate scaling factor is incorporated
into the result. For example, ``Hz`` and ``s`` will now appropriately be recognized
as inverses:

    >>> from yt.units import Hz, s
    >>> frequency = 60*Hz
    >>> time = 60*s
    >>> frequency*time
    unyt_quantity(3600, '(dimensionless)')

Similar simplifications will happen even if units aren't reciprocals of each
other, for example here ``hour`` and ``minute`` automatically cancel each other:

    >>> from yt.units import erg, minute, hour
    >>> power = [20, 40, 80] * erg / minute
    >>> elapsed_time = 3*hour
    >>> print(power*elapsed_time)
    [ 3600.  7200. 14400.] erg

Alternate Unit Name Resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It's now possible to use a number of common alternate spellings for unit names
and if ``unyt`` knows about the alternate spelling it will automatically resolve
alternate spellings to a canonical name. For example, it's now possible to do
things like this:

    >>> import yt.units as u
    >>> d = 20*u.mile
    >>> d.to('km')
    unyt_quantity(32.18688, 'km')
    >>> d.to('kilometer')
    unyt_quantity(32.18688, 'km')
    >>> d.to('kilometre')
    unyt_quantity(32.18688, 'km')

You can also use alternate unit names in more complex algebraic unit expressions:

    >>> v = d / (20*u.minute)
    >>> v.to('kilometre/hour')
    unyt_quantity(96.56064, 'km/hr')

In this example the common british spelling ``"kilometre"`` is resolved to
``"km"`` and ``"hour"`` is resolved to ``"hr"``.

New Method for Accessing Sample Datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is now a function entitled ``load_sample()`` that allows the user to
automatically load sample data from the yt hub in a local yt session.
Previously, users would have to explicitly download these data directly from
`https://yt-project.org/data <https://yt-project.org/data>`_, unpackage them,
and load them into a yt session, but now this occurs from within a python
session.  For more information see:
:ref:`Loading Sample Data <loading-sample-data>`

API Changes
-----------

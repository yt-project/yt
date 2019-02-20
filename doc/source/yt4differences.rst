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

Background
----------

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

Updating to yt 4.0 from Old Versions (and going back)
-----------------------------------------------------


.. _transitioning-to-4.0:

Converting Old Scripts to Work with yt 4.0
------------------------------------------


Cool New Things
---------------


Scatter and gather approach for SPH data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

    ds = yt.load('snapshot_033/snap_033.0.hdf5')

    plot = yt.SlicePlot(ds, 2, ('gas', 'density'))
    plot.save()

    plot = yt.ProjectionPlot(ds, 2, ('gas', 'density'))
    plot.save()

    arbitrary_grid = ds.arbitrary_grid([0.0, 0.0, 0.0], [25, 25, 25],
                                       dims=[16, 16, 16])
    ag_density = arbitrary_grid[('gas', 'density')]

    covering_grid = ds.covering_grid(4, 0, 16)
    cg_density = covering_grid[('gas', 'density')]

In the above example the `covering_grid` and the `arbitrary_grid` will return
the same data. In fact, these containers are very similar but provide a
slighlty different API.

The above code can be modified to use the gather approach by changing a global
setting for the dataset. This can be achieved with
`ds.sph_smoothing_style = "gather"`, so far, the gather approach is not
supported for projections.

The default behaviour for SPH interpolation is that the values are normalized
inline with Eq. 9 in `SPLASH, Price (2009) <https://arxiv.org/pdf/0709.0832.pdf>`_.
This can be disabled with `ds.use_sph_normalization = False`. This will
disable the normalization for all future interpolations.

The gather approach requires finding nearest neighbors using the KDTree. The
first call will generate a KDTree for the entire dataset which will be stored in
a sidecar file. This will be loaded whenever neccesary.

Off-Axis Projection for SPH Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The current `OffAxisProjectionPlot` class will now support SPH projection plots.

The following is a code example:

.. code-block:: python

    import yt

    ds = yt.load('Data/GadgetDiskGalaxy/snapshot_200.hdf5')

    smoothing_field = ('gas', 'density')

    _, center = ds.find_max(smoothing_field)

    sp = ds.sphere(center, (10, 'kpc'))

    normal_vector = sp.quantities.angular_momentum_vector()

    prj = yt.OffAxisProjectionPlot(ds, normal_vector, smoothing_field, center, (20, 'kpc'))

    prj.save()

Octree
^^^^^^

Whilst the move away from the global octree is a promising one in terms of
perfomance and dealing with SPH data in a more intuitive manner, it does remove
a useful feature. We are aware that many uses will have older scripts which take
advantage of the global octree.

As such, we have added support to smooth SPH data onto an octree when desired by
the users. The new octree is designed to give results consistent with those of
the previous octree, but the new octree takes advantage of the scatter and
gather machinery also added.

It should be noted that the 

.. code-block:: python

    import yt

    ds = yt.load('GadgetDiskGalaxy/snapshot_200.hdf5')
    left = np.array([0, 0, 0], dtype='float64')
    right = np.array([64000, 64000, 64000], dtype='float64')

    # generate an octree
    octree = ds.octree(left, right, n_ref=64)

    # the density will be calculated using SPH scatter
    density = octree[('PartType0', 'density')]

    # this will return the x positions of the octs
    x = octree[('index', 'x')]

The above code can be modified to use the scatter approach by using
`ds.sph_smoothing_style = 'gather'` before any field access. The octree also
accepts `over_refine_factor` which behaves identically to that in the previous
branch, describing how many cells are in each leaf.

The new octree also has the ability to not be an octree. We have a new kwarg,
`density_factor` which allows the construction of dense trees. In a traditional
octree, if a leaf has more particles that a critical value `n_ref`, then it
divides into 8 new children (hence the name oct). The value of `density_factor`
allows the node to divide into 2^(3*density_factor).

API Changes
-----------

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

It is now possible to generate slice plots, projection plots and arbitrary grids
of smoothed quanitities using these operations. The following code demonstrates
how this could be achieved:

.. code-block:: python

    import yt

    ds = yt.load('snapshot_033/snap_033.0.hdf5')

    plot = yt.SlicePlot(ds, 2, ('gas', 'density'))
    plot.save()

    plot = yt.ProjectionPlot(ds, 2, ('gas', 'density'))
    plot.save()

    arbitrary_grid = ds.arbitrary_grid([0.0, 0.0, 0.0], [5, 5, 5],
                                       dims=[10, 10, 10])
    density = arbitrary_grid[('gas', 'density')]

The default behaviour for sPH interpolation is that the values are normalized
inline with Eq. 9 in `SPLASH, Price (2009) <https://arxiv.org/pdf/0709.0832.pdf>`_.
This can be disabled with `ds.use_sph_normalization = False`.

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


API Changes
-----------


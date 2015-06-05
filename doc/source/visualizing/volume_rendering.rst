.. _volume_rendering:

3D Visualization and Volume Rendering
=====================================

yt has the ability to create 3D visualizations, using a process known as volume
rendering.  At present, its capabilities are limited to software volume
rendering, but porting this capability to OpenGL and GPUs is an area of active
development.

Constructing a 3D visualization is a process of describing the "scene" that
will be rendered.  This includes the location of the viewing point (i.e., where
the "camera" is placed), the method by which a system would be viewed (i.e.,
the "lens," which may be orthographic, perspective, fisheye, spherical, and so
on) and the components that will be rendered (render "sources," such as volume
elements, lines, annotations, and opaque surfaces).  The 3D plotting
infrastructure then develops a resultant image from this scene, which can be
saved to a file or viewed inline.

By constructing the scene in this programmatic way, full control can be had
over each component in the scene as well as the method by which the scene is
rendered; this can be used to prototype visualizations, inject annotation such
as grid or continent lines, and then to render a production-quality
visualization.  By changing the "lens" used, a single camera path can output
images suitable for planetarium domes, immersive and head tracking systems
(such as the Occulus Rift or recent "spherical" movie viewers such as the
mobile YouTube app), as well as standard screens.

.. image:: _images/scene_diagram.svg
   :width: 50%
   :align: center
   :alt: Diagram of a 3D Scene

In versions of yt prior to 3.2, the only volume rendering interface accessible
was through the "camera" object.  This presented a number of problems,
principle of which was the inability to describe new scene elements or to
develop complex visualizations that were independent of the specific elements
being rendered.  The new "scene" based interface present in yt 3.2 and beyond
enables both more complex visualizations to be constructed as well as a new,
more intuitive interface for very simple 3D visualizations.

.. warning:: 3D visualizations can be fun but frustrating!  Tuning the
             parameters to both look nice and convey useful scientific
             information can be hard.  We've provided information about best
             practices and tried to make the interface easy to develop nice
             visualizations, but getting them *just right* is often
             time-consuming.

Tutorial
--------

The scene interface provides a more modular interface for creating renderings
of arbitrary data sources. As such, manual composition of a scene can require a
bit more work, but we will also provide several helper functions that attempt
to create satisfactory default volume renderings.

.. note:: It's usually best to start out simple with the built-in helper
          interface, and expand on that if you need to.

Here is a working example for rendering the IsolatedGalaxy dataset.

.. python-script::
  import yt
  # load the data
  ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
  # volume render the 'density' field, and save the resulting image
  im, sc = yt.volume_rendering(ds, 'density', fname='test_rendering.png')

  # im is the image that was generated.
  # sc is an instance of a Scene object, which allows you to further refine
  # your renderings.

When the volume_rendering function is called, first an empty
:class:`~yt.visualization.volume_rendering.scene.Scene` object is
created. Next, a 
:class:`~yt.visualization.volume_rendering.api.VolumeSource`
object is created, which deomposes the volume elements
into an tree structure to provide back-to-front rendering of fixed-resolution
blocks of data.  (If the volume elements are grids, this uses a
:class:`~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree` object.) When the
:class:`~yt.visualization.volume_rendering.api.VolumeSource`
object is created, by default it will create a transfer function
based on the extrema of the field that you are rendering. The transfer function
describes how rays that pass through the domain are "transfered" and thus how
brightness and color correlates to the field values.  Modifying and adjusting
the transfer function is the primary way to modify the appearance of an image
based on volumes.

Once the basic set of objects to be rendered is constructed, a
:class:`~yt.visualization.volume_rendering.camera.Camera` object is created and
added to the scene.  By default the creation of a camera also creates a
default, plane-parallel :class:`~yt.visualization.volume_rendering.lens.Lens`
object. The analog to a real camera is intentional -- a camera can take a
picture of a scene from a particular point in time and space, but different
lenses can be swapped in and out.  For example, this might include a fisheye
lens, a spherical lens, or some other method of describing the direction and
origin of rays for rendering. Once the camera is added to the scene object, we
call the main method of the
:class:`~yt.visualization.volume_rendering.scene.Scene` class,
:meth:`~yt.visualization.volume_rendering.scene.Scene.render` .  When called,
the scene will loop through all of the
:class:`~yt.visualization.volume_rendering.render_source.RenderSource` objects
that have been added, and integrate the radiative transfer equation through the
volume. Finally, the image and scene object is returned to the user.

In this example, we don't add on any non-volume rendering sources; however, if
such sources are added, they will be integrated as well.

Modifying the Scene
-------------------

Once a basic scene has been created, where default render sources and basic
camera operations are utilized, deeper modifications are possible.  These
modifications can tune the appearance of the render sources (such as which
colors correspond to which values in the data) as well as the shape of the
rendered image, the position of the camera in the scene, and other elements
present in the scene.  Below, we describe a few of the aspects of tuning a
scene to create a visualization that is communicative and pleasing.

.. _transfer_functions:

Transfer Functions
++++++++++++++++++

Transfer functions are the most essential component of a rendering that
includes volume sources.  Several different fundamental types have been
provided, but there are many different ways the construct complicated
expressions to produce visualizations and images using the underlying
machinery.

.. note::
   All of the information about how transfer functions are used and values
   extracted is contained in the functions `TransferFunctionProxy.eval_transfer`
   and `FIT_get_value` in the file `yt/_amr_utils/VolumeIntegrator.pyx`.  If
   you're curious about how to construct your own, or why you get the values
   you do, you should read the source!

There are three ready-to-go transfer functions implemented in yt.
:class:`~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction`,
:class:`~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction`,
and
:class:`~yt.visualization.volume_rendering.transfer_functions.PlanckTransferFunction`.

Color Transfer Functions
^^^^^^^^^^^^^^^^^^^^^^^^

These transfer functions are the standard way to apply colors to specific
values in the field being rendered.  For instance, applying isocontours at
specific densities.  They have several different mechanisms that can be used.
The easiest mechanism is to use
:meth:`~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction.add_layers`,
which will add evenly spaced isocontours between the bounds of the transfer
function.  However, you can also use
:meth:`~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction.sample_colormap`,
which will sample a colormap at a given value.  Additionally, you can directly
call
:meth:`~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction.add_gaussian`,
which will allow you to specify the colors directly.

An alternate method for modifying the colormap is done using
:meth:`~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction.map_to_colormap`,
where you can map a segment of the transfer function space to an entire
colormap at a single alpha value.  This is sometimes useful for very opaque
renderings.

See :ref:`cookbook-simple_volume_rendering` for an example usage.

Projection Transfer Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is designed to allow you to very easily project off-axis through a region.
See :ref:`cookbook-offaxis_projection` for a simple example.  Note that the
integration here is scaled to a width of 1.0; this means that if you want to
apply a colorbar, you will have to multiply by the integration width (specified
when you initialize the volume renderer) in whatever units are appropriate.

Planck Transfer Function
^^^^^^^^^^^^^^^^^^^^^^^^

This transfer function is designed to apply a semi-realistic color field based
on temperature, emission weighted by density, and approximate scattering based
on the density.  This class is currently under-documented, and it may be best
to examine the source code to use it.

More Complicated Transfer Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For more complicated transfer functions, you can use the
:class:`~yt.visualization.volume_rendering.transfer_functions.MultiVariateTransferFunction`
object.  This allows for a set of weightings, linkages and so on.

.. _transfer-function-helper:

TransferFunctionHelper
----------------------

.. notebook:: TransferFunctionHelper_Tutorial.ipynb

Adding New Sources
++++++++++++++++++

The resulting image of a rendering process is a combination of the different
sources present in a scene.  While at present there are only a few sources
available, in principle new sources can be defined and added to yt over time.

By default, the scene will construct a volume object that includes the fluid
components of a data source. 

Volume Objects
++++++++++++++

Hard and Opaque Objects
+++++++++++++++++++++++

Annotations
+++++++++++

Moving Around
-------------

.. _camera_movement:

Moving the Camera
+++++++++++++++++

There are multiple ways to manipulate the camera viewpoint to create a series of
renderings.  For an example, see this cookbook:
:ref:`cookbook-camera_movement`.  For a current list of
motion helper functions, see the docstrings associated with
:class:`~yt.visualization.volume_rendering.camera.Camera`.

Looking at Things
+++++++++++++++++

Orienting the Camera
++++++++++++++++++++

Changing Lenses
+++++++++++++++

Volume Rendering Method
-----------------------

Direct ray casting through a volume enables the generation of new types of
visualizations and images describing a simulation.  yt has the facility
to generate volume renderings by a direct ray casting method.  However, the
ability to create volume renderings informed by analysis by other mechanisms --
for instance, halo location, angular momentum, spectral energy distributions --
is useful.

The volume rendering in yt follows a relatively straightforward approach.

#. Create a set of transfer functions governing the emission and absorption as
   a function of one or more variables. (:math:`f(v) \rightarrow (r,g,b,a)`)
   These can be functions of any field variable, weighted by independent
   fields, and even weighted by other evaluated transfer functions.  (See
   `transfer_functions`.)
#. Partition all chunks into non-overlapping, fully domain-tiling "bricks."
   Each of these "bricks" contains the finest available data at any location.
#. Generate vertex-centered data for all grids in the volume rendered domain.
#. Order the bricks from front-to-back.
#. Construct plane of rays parallel to the image plane, with initial values set
   to zero and located at the back of the region to be rendered.
#. For every brick, identify which rays intersect.  These are then each 'cast'
   through the brick.

   #. Every cell a ray intersects is sampled 5 times (adjustable by parameter),
      and data values at each sampling point are trilinearly interpolated from
      the vertex-centered data.
   #. Each transfer function is evaluated at each sample point.  This gives us,
      for each channel, both emission (:math:`j`) and absorption
      (:math:`\alpha`) values.
   #. The value for the pixel corresponding to the current ray is updated with
      new values calculated by rectangular integration over the path length:

      :math:`v^{n+1}_{i} =  j_{i}\Delta s + (1 - \alpha_{i}\Delta s )v^{n}_{i}`

      where :math:`n` and :math:`n+1` represent the pixel before and after
      passing through a sample, :math:`i` is the color (red, green, blue) and 
      :math:`\Delta s` is the path length between samples.
   #. Determine if any addition integrate will change the sample value; if not,
      terminate integration.  (This reduces integration time when rendering
      front-to-back.)
#. The image is returned to the user:

.. image:: _images/vr_sample.jpg
   :width: 512

MPI Parallelization
-------------------
Currently the volume renderer is parallelized using MPI to decompose the volume
by attempting to split up the
:class:`~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree` in a balanced way.  This
has two advantages: 

#.  The :class:`~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree`
    construction is parallelized since each MPI task only needs
    to know about the part of the tree it will traverse.
#.  Each MPI task will only read data for portion of the volume that it has
    assigned.

Once the :class:`~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree` has been 
constructed, each MPI task begins the rendering
phase until all of its bricks are completed.  At that point, each MPI task has
a full image plane which we then use a tree reduction to construct the final
image, using alpha blending to add the images together at each reduction phase.

Caveats:

#.  At this time, the :class:`~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree`
    can only be decomposed by a power of 2 MPI
    tasks.  If a number of tasks not equal to a power of 2 are used, the largest
    power of 2 below that number is used, and the remaining cores will be idle.
    This issue is being actively addressed by current development.
#.  Each MPI task, currently, holds the entire image plane.  Therefore when
    image plane sizes get large (>2048^2), the memory usage can also get large,
    limiting the number of MPI tasks you can use.  This is also being addressed
    in current development by using image plane decomposition.

For more information about enabling parallelism, see :ref:`parallel-computation`.

OpenMP Parallelization
----------------------

The volume rendering also parallelized using the OpenMP interface in Cython.
While the MPI parallelization is done using domain decomposition, the OpenMP
threading parallelizes the rays intersecting a given brick of data.  As the
average brick size relative to the image plane increases, the parallel
efficiency increases. 

By default, the volume renderer will use the total number of cores available on
the symmetric multiprocessing (SMP) compute platform.  For example, if you have
a shiny new laptop with 8 cores, you'll by default launch 8 OpenMP threads.
The number of threads can be controlled with the num_threads keyword in
:meth:`~yt.visualization.volume_rendering.camera.Camera.snapshot`.  You may also restrict the number of OpenMP threads used
by default by modifying the environment variable OMP_NUM_THREADS. 

Running in Hybrid MPI + OpenMP
------------------------------

The two methods for volume rendering parallelization can be used together to
leverage large supercomputing resources.  When choosing how to balance the
number of MPI tasks vs OpenMP threads, there are a few things to keep in mind.
For these examples, we will assume you are using Nmpi MPI tasks, and Nmp OpenMP
tasks, on a total of P cores. We will assume that the machine has a Nnode SMP
nodes, each with cores_per_node cores per node.

#.  For each MPI task, num_threads (or OMP_NUM_THREADS) OpenMP threads will be
    used. Therefore you should usually make sure that Nmpi*Nmp = P.  
#.  For simulations with many grids/AMRKDTree bricks, you generally want to increase Nmpi.
#.  For simulations with large image planes (>2048^2), you generally want to
    decrease Nmpi and increase Nmp. This is because, currently, each MPI task
    stores the entire image plane, and doing so can approach the memory limits
    of a given SMP node.
#.  Please make sure you understand the (super)computer topology in terms of
    the numbers of cores per socket, node, etc when making these decisions.
#.  For many cases when rendering using your laptop/desktop, OpenMP will
    provide a good enough speedup by default that it is preferable to launching
    the MPI tasks.

For more information about enabling parallelism, see :ref:`parallel-computation`.

Opacity
-------

There are currently two models for opacity when rendering a volume, which are
controlled in the ColorTransferFunction with the keyword
grey_opacity=False(default)/True. The first (default) will act such for each of
the r,g,b channels, each channel is only opaque to itself.  This means that if
a ray that has some amount of red then encounters material that emits blue, the
red will still exist and in the end that pixel will be a combination of blue
and red.  However, if the ColorTransferFunction is set up with
grey_opacity=True, then blue will be opaque to red, and only the blue emission
will remain.  

For an in-depth example, please see the cookbook example on opaque renders here: 
:ref:`cookbook-opaque_rendering`.

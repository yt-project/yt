.. _volume_rendering:

Volume Rendering: Making 3D Photorealistic Isocontoured Images
==============================================================

Volume rendering, as implemented in yt, is a mechanism by which rays are cast
through a domain, converting field values to emission and absorption, and producing a final image.
This provides the ability to create off-axis projections, isocontour images,
volume emission, and absorption from intervening material.  The primary goal 
of the volume rendering in yt is to provide the ability to make
*scientifically-informed* visualizations of simulations.  

The volume renderer is implemented in a hybrid of Python and Cython, which is
Python-like code compiled down to C.  It has been optimized, but it is still a
*software* volume renderer: it does not currently operate on graphics
processing units (GPUs).  However, while the rendering engine itself may not
directly translate to GPU code (OpenCL, CUDA or OpenGL), the Python structures:
partitioning, transfer functions, display, etc., may be useful in the future
for transitioning the rendering to the GPU.  In addition, this allows users to create
volume renderings on traditional supercomputing platforms that may not have access to GPUs.

The volume renderer is also threaded using OpenMP.  Many of the commands
(including `snapshot`) will accept a `num_threads` option.

Tutorial
--------

Volume renderings are created by combining three objects: a volume
homogenization; a transfer function, and a camera object.

#. Find the appropriate bounds for your data.
#. Create a ColorTransferFunction object.
#. Create a Camera object, which homogenizes the volume and orients the viewing
   direction
#. Take a snapshot and save the image.

Here is a working example for the IsolatedGalaxy dataset.

.. python-script::

   import yt
   import numpy as np

   ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
   # Choose a field
   field = 'density'
   # Do you want the log of the field?
   use_log = True

   # Find the bounds in log space of for your field
   dd = ds.all_data()
   mi, ma = dd.quantities.extrema(field)

   if use_log:
       mi,ma = np.log10(mi), np.log10(ma)

   # Instantiate the ColorTransferfunction.
   tf = yt.ColorTransferFunction((mi, ma))

   # Set up the camera parameters: center, looking direction, width, resolution
   c = (ds.domain_right_edge + ds.domain_left_edge)/2.0
   L = np.array([1.0, 1.0, 1.0])
   W = ds.quan(0.3, 'unitary')
   N = 256 

   # Create a camera object
   cam = ds.camera(c, L, W, N, tf, fields = [field], log_fields = [use_log])

   # Now let's add some isocontours, and take a snapshot, saving the image
   # to a file.
   tf.add_layers(10, 0.01, colormap = 'RdBu_r')
   im = cam.snapshot('test_rendering.png')

   # To add the domain box to the image:
   nim = cam.draw_domain(im)
   nim.write_png('test_rendering_with_domain.png')

   # To add the grid outlines to the image:
   nim = cam.draw_grids(im)
   nim.write_png('test_rendering_with_grids.png')

Method
------

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
#. Order the bricks from back-to-front.
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
#. The image is returned to the user:

.. image:: _images/vr_sample.jpg
   :width: 512

.. _the-camera-interface:

The Camera Interface
--------------------

A camera object has also been created, to allow for more programmatic
descriptions of the viewpoint and image plane, and to allow for moving the
camera object through the volume and creating multiple images.  There are
several camera objects available, but the most commonly used is the standard,
orthographic projection camera.

The primary interface here is through the creation of an instance of
:class:`~yt.visualization.volume_rendering.camera.Camera`, which represents a
viewpoint into a volume.  The camera optionally accepts a volume, which can be
either an instance of
:class:`~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree` that
has already been initialized.  If one is not supplied, the camera will generate
one itself.  This can also be specified if you wish to save bricks between
repeated calls, thus saving considerable amounts of time.

The camera interface allows the user to move the camera about the domain, as
well as providing interfaces for zooming in and out.  Furthermore, yt now
includes a stereoscopic camera
(:class:`~yt.visualization.volume_rendering.camera.StereoPairCamera`).

Much like most data objects, the
:class:`~yt.visualization.volume_rendering.camera.Camera` object hangs off of
the index file, and can be instantiated in that manner.

.. warning::  The keyword *no_ghost* has been set to True by default
              for speed considerations.  However, because this turns off ghost
              zones, there may be artifacts at grid boundaries.  If a higher quality
              rendering is required, use *no_ghost = False*.

Here's a fully functional script that demonstrates how to use the camera
interface.

For an example, see the cookbook :ref:`cookbook-simple_volume_rendering`.

The :class:`~yt.visualization.volume_rendering.camera.StereoPairCamera` object
has a single primary method,
:meth:`~yt.visualization.volume_rendering.camera.StereoPairCamera.split`, that
will return two cameras, a left and a right.

.. _camera_movement:

Camera Movement
---------------

There are multiple ways to manipulate the camera viewpoint to create a series of
renderings.  For an example, see this cookbook:
:ref:`cookbook-camera_movement`.  For a current list of
motion helper functions, see the docstrings associated with
:class:`~yt.visualization.volume_rendering.camera.Camera`.

.. _transfer_functions:

Transfer Functions
------------------

Transfer functions are the most essential component.  Several different
fundamental types have been provided, but there are many different ways the
construct complicated expressions to produce visualizations and images using
the underlying machinery.

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

.. _healpix_volume_rendering:

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

Lighting
--------

Lighting can be optionally used in volume renders by specifying use_light=True
in the Camera object creation.  If used, one can then change the default
lighting color and direction by modifying Camera.light_dir and
Camera.light_rgb.  Lighting works in this context by evaluating not only the
field value but also its gradient in order to compute the emissivity.  This is
not the same as casting shadows, but provides a way of highlighting sides of a
contour.  

Generating a Homogenized Volume
-------------------------------

In order to perform a volume rendering, the data must first be decomposed into
a HomogenizedVolume object.  This structure splits the domain up into
single-resolution tiles which cover the domain at the highest resolution
possible for a given point in space.  This means that every point in space is
mapped to exactly one data point, which receives its values from the highest
resolution grid that covers that volume.

The creation of these homogenized volumes is done during the 
:class:`~yt.visualization.volume_rendering.camera.Camera`  object
instantiation by default.  However, in some cases it is useful to first build
your homogenized volume to then be passed in to the camera. A sample usage is shown
in :ref:`cookbook-amrkdtree_downsampling`.

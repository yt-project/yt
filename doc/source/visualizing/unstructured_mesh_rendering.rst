.. _unstructured_mesh_rendering:

Unstructured Mesh Rendering
===========================

Beginning with version 3.3, yt has the ability to volume render unstructured
meshes from, for example, finite element calculations. In order to use this
capability, a few additional dependencies are required beyond those you get
when you run the install script. First, `embree <https://embree.github.io>`
(a fast software ray-tracing library from Intel) must be installed, either
by compiling from source or by using one of the pre-built binaries available
at Embree's `downloads <https://embree.github.io/downloads.html>` page. Once
Embree is installed, you must also create a symlink next to the library. For
example, if the libraries were installed at /usr/local/lib/, you must do

.. code-block:: bash

    sudo ln -s /usr/local/lib/libembree.2.6.1.dylib /usr/local/lib/libembree.so

Second, the python bindings for embree (called 
`pyembree <https://github.com/scopatz/pyembree>`) must also be installed. To
do so, first obtain a copy, by .e.g. cloning the repo:

.. code-block:: bash

    git clone https://github.com/scopatz/pyembree

To install, navigate to the root directory and run the setup script:

.. code-block:: bash

    python setup.py develop

If Embree was installed to some location that is not in your path by default,
you will need to pass in CFLAGS and LDFLAGS to the setup.py script. For example,
the Mac OS package installer puts the installation at /opt/local/ instead of 
usr/local. To account for this, you would do:

.. code-block:: bash

    CFLAGS='-I/opt/local/include' LDFLAGS='-L/opt/local/lib' python setup.py install

You must also use these flags when building any part of yt that links against
pyembree.

Once the pre-requisites are installed, unstructured mesh data can be rendered
much like any other dataset. In particular, a new type of 
:class:`~yt.visualization.volume_rendering.render_source.RenderSource` object
has been defined, called the 
:class:`~yt.visualization.volume_rendering.render_source.MeshSource`, that
represents the unstructured mesh data that will be rendered. The user creates 
this object, and also defines a
:class:`~yt.visualization.volume_rendering.camera.Camera` 
that specifies your viewpoint into the scene. When 
:class:`~yt.visualization.volume_rendering.render_source.RenderSource` is called,
a set of rays are cast at the source. Each time a ray strikes the source mesh,
the data is sampled at the intersection point at the resulting value gets 
saved into an image.

See below for examples. First, here is an example of rendering a hexahedral mesh.

.. python-script::
   import yt
   import pylab as plt
   from yt.visualization.volume_rendering.render_source import MeshSource
   from yt.visualization.volume_rendering.camera import Camera
   from yt.utilities.exodusII_reader import get_data

   # load the data
   coords, connectivity, data = get_data("data/out.e-s010")
   mesh_id = 0
   field_name = ('gas', 'diffused')
   ds = yt.load_unstructured_mesh(data[mesh_id], connectivity[mesh_id], coords[mesh_id])

   # create the RenderSource
   ms = MeshSource(ds, field_name)

   # set up camera
   cam = Camera(ds)
   camera_position = ds.arr([-3.0, 3.0, -3.0], 'code_length')
   north_vector = ds.arr([0.0, 1.0, 0.0], 'dimensionless')
   cam.resolution = (800, 800)
   cam.set_position(camera_position, north_vector)

   # make the image
   im = ms.render(cam)

   # plot and save
   plt.imshow(im, cmap='Eos A', origin='lower', vmin=0, vmax=2.0)
   plt.gca().axes.get_xaxis().set_visible(False)
   plt.gca().axes.get_yaxis().set_visible(False)
   cb = plt.colorbar()
   cb.set_label(field_name[1])
   plt.savefig('hex_mesh_render.png')

Next, here is an example of rendering a dataset with tetrahedral mesh elements.

.. python-script::
   import yt
   import pylab as plt
   from yt.visualization.volume_rendering.render_source import MeshSource
   from yt.visualization.volume_rendering.camera import Camera
   from yt.utilities.exodusII_reader import get_data

   # load the data
   filename = "../moose/test/tests/mesh/high_order_elems/gold/high_order_elems_tet4_refine_out.e"
   coords, connectivity, data = get_data(filename)
   mesh_id = 0
   field_name = ('gas', 'u')
   ds = yt.load_unstructured_mesh(data[mesh_id], connectivity[mesh_id], coords[mesh_id])

   # create the RenderSource
   ms = MeshSource(ds, field_name)

   # set up camera
   cam = Camera(ds)
   camera_position = ds.arr([3.0, 3.0, 3.0], 'code_length')
   cam.set_width(ds.arr([2.0, 2.0, 2.0], 'code_length'))
   north_vector = ds.arr([0.0, 1.0, 0.0], 'dimensionless')
   cam.resolution = (800, 800)
   cam.set_position(camera_position, north_vector)

   # make the image
   im = ms.render(cam)

   # plot and save
   plt.imshow(im, cmap='Eos A', origin='lower', vmin=0.0, vmax=1.0)
   plt.gca().axes.get_xaxis().set_visible(False)
   plt.gca().axes.get_yaxis().set_visible(False)
   cb = plt.colorbar()
   cb.set_label(field_name[1])
   plt.savefig('tet_mesh_render.png')

Finally, here is a script that creates frames of a movie. It calls the rotate()
method 300 times, saving a new image to the disk each time.

.. python-script::
   import yt
   import pylab as plt
   from yt.visualization.volume_rendering.render_source import MeshSource
   from yt.visualization.volume_rendering.camera import Camera
   from yt.utilities.exodusII_reader import get_data

   # load dataset
   coords, connectivity, data = get_data("data/out.e-s010")
   mesh_id = 0
   field_name = ('gas', 'diffused')
   ds = yt.load_unstructured_mesh(data[mesh_id], connectivity[mesh_id], coords[mesh_id])

   # create the RenderSource
   ms = MeshSource(ds, field_name)

   # set up camera
   cam = Camera(ds)
   camera_position = ds.arr([-3.0, 3.0, -3.0], 'code_length')
   north_vector = ds.arr([0.0, 1.0, 0.0], 'dimensionless')
   cam.set_position(camera_position, north_vector)
   cam.steady_north = True

   # make movie frames
   num_frames = 301
   for i in range(num_frames):
       cam.rotate(2.0*np.pi/num_frames)
       im = ms.render(cam)
       plt.imshow(im, cmap='Eos A', origin='lower',vmin=0.0, vmax=2.0)
       plt.gca().axes.get_xaxis().set_visible(False)
       plt.gca().axes.get_yaxis().set_visible(False)
       cb = plt.colorbar()
       cb.set_label('diffused')
       plt.savefig('movie_frames/surface_render_%.4d.png' % i)
       plt.clf()

.. _unstructured_mesh_rendering:

Unstructured Mesh Rendering
===========================

Beginning with version 3.3, yt has the ability to volume render unstructured
meshes from, for example, finite element calculations. In order to use this
capability, a few additional dependencies are required beyond those you get
when you run the install script. First, embree (a fast software ray-tracing
library from Intel) must be installed, following the instructions here. 
Second, the python bindings for embree (called ''pyembree'') must also 
be installed. 

Once the pre-requisites are installed, unstructured mesh data can be rendered
much like any other dataset. In particular, a new type of RenderSource object
has been defined, called the MeshSource, that represents the unstructured mesh
data that will be rendered. The user creates this object, and also defines a 
camera that specifies your viewpoint into the scene. When render() is called,
a set of rays are cast at the source. Each time a ray strikes the source mesh,
the data is sampled at the intersection point at the resulting value gets 
saved into an image.

See below for examples. First, here is an example of rendering a hexahedral mesh.

Next, here is an example of rendering a dataset with tetrahedral mesh elements.

Finally, here is a script that creates frames of a movie. It calls the rotate()
method 300 times, saving a new image to the disk each time.

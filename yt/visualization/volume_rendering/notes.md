
Overview of Volume Rendering
============================

In 3.0, we have moved away from the "god class" that was Camera, and have
attempted to break down the VR system into a hierarchy of classes.  So far
we are at:

1. Scene 
2. Camera
3. Lens 
4. Source

For now, a scene only has one camera, i.e. one viewpoint. I would like this to be
extended to multiple cameras at some point, but not in this pass.

A Camera can have many lenses. When taking a snapshot, the Camera will loop 
over the lenses that have been added by the user.  We should come up with a
naming convention and storage system.


A Lens defines how the vectors are oriented pointing outward from the camera
position.  Plane-parallel, Perspective, Fisheye are the first set that need to
be implemented. As much of the Lens as possible will be set up using defaults 
derived from the scene, such as the width/depth/etc.

A Source is a data source with intent on how to visualize it.  For example, a
VolumeSource should be treated volumetrically, with a transfer function defined
for a given field or set of fields.  A generic OpaqueSource should define
a method for pixelizing a ZBuffer object, carrying information about both the
color and depth of the surface/streamline/annotation. These will be used for
compositing later.


sc = Scene(data_source)
cam = sc.add_camera(cam) // triggers cam.set_defaults_from_data_source(data_source)
lens = PlaneParallelLens()
cam.set_lens(lens) # This sets up lens based on camera.


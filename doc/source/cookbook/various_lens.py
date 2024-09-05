import numpy as np

import yt
from yt.visualization.volume_rendering.api import Scene, create_volume_source

field = ("gas", "density")

# normal_vector points from camera to the center of the final projection.
# Now we look at the positive x direction.
normal_vector = [1.0, 0.0, 0.0]
# north_vector defines the "top" direction of the projection, which is
# positive z direction here.
north_vector = [0.0, 0.0, 1.0]

# Follow the simple_volume_rendering cookbook for the first part of this.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
sc = Scene()
vol = create_volume_source(ds, field=field)
tf = vol.transfer_function
tf.grey_opacity = True

# Plane-parallel lens
cam = sc.add_camera(ds, lens_type="plane-parallel")
# Set the resolution of the final projection.
cam.resolution = [250, 250]
# Set the location of the camera to be (x=0.2, y=0.5, z=0.5)
# For plane-parallel lens, the location info along the normal_vector (here
# is x=0.2) is ignored.
cam.position = ds.arr(np.array([0.2, 0.5, 0.5]), "code_length")
# Set the orientation of the camera.
cam.switch_orientation(normal_vector=normal_vector, north_vector=north_vector)
# Set the width of the camera, where width[0] and width[1] specify the length and
# height of final projection, while width[2] in plane-parallel lens is not used.
cam.set_width(ds.domain_width * 0.5)
sc.add_source(vol)
sc.save("lens_plane-parallel.png", sigma_clip=6.0)

# Perspective lens
cam = sc.add_camera(ds, lens_type="perspective")
cam.resolution = [250, 250]
# Standing at (x=0.2, y=0.5, z=0.5), we look at the area of x>0.2 (with some open angle
# specified by camera width) along the positive x direction.
cam.position = ds.arr([0.2, 0.5, 0.5], "code_length")
cam.switch_orientation(normal_vector=normal_vector, north_vector=north_vector)
# Set the width of the camera, where width[0] and width[1] specify the length and
# height of the final projection, while width[2] specifies the distance between the
# camera and the final image.
cam.set_width(ds.domain_width * 0.5)
sc.add_source(vol)
sc.save("lens_perspective.png", sigma_clip=6.0)

# Stereo-perspective lens
cam = sc.add_camera(ds, lens_type="stereo-perspective")
# Set the size ratio of the final projection to be 2:1, since stereo-perspective lens
# will generate the final image with both left-eye and right-eye ones jointed together.
cam.resolution = [500, 250]
cam.position = ds.arr([0.2, 0.5, 0.5], "code_length")
cam.switch_orientation(normal_vector=normal_vector, north_vector=north_vector)
cam.set_width(ds.domain_width * 0.5)
# Set the distance between left-eye and right-eye.
cam.lens.disparity = ds.domain_width[0] * 1.0e-3
sc.add_source(vol)
sc.save("lens_stereo-perspective.png", sigma_clip=6.0)

# Fisheye lens
dd = ds.sphere(ds.domain_center, ds.domain_width[0] / 10)
cam = sc.add_camera(dd, lens_type="fisheye")
cam.resolution = [250, 250]
v, c = ds.find_max(field)
cam.set_position(c - 0.0005 * ds.domain_width)
cam.switch_orientation(normal_vector=normal_vector, north_vector=north_vector)
cam.set_width(ds.domain_width)
cam.lens.fov = 360.0
sc.add_source(vol)
sc.save("lens_fisheye.png", sigma_clip=6.0)

# Spherical lens
cam = sc.add_camera(ds, lens_type="spherical")
# Set the size ratio of the final projection to be 2:1, since spherical lens
# will generate the final image with length of 2*pi and height of pi.
# Recommended resolution for YouTube 360-degree videos is [3840, 2160]
cam.resolution = [500, 250]
# Standing at (x=0.4, y=0.5, z=0.5), we look in all the radial directions
# from this point in spherical coordinate.
cam.position = ds.arr([0.4, 0.5, 0.5], "code_length")
cam.switch_orientation(normal_vector=normal_vector, north_vector=north_vector)
# In (stereo)spherical camera, camera width is not used since the entire volume
# will be rendered
sc.add_source(vol)
sc.save("lens_spherical.png", sigma_clip=6.0)

# Stereo-spherical lens
cam = sc.add_camera(ds, lens_type="stereo-spherical")
# Set the size ratio of the final projection to be 1:1, since spherical-perspective lens
# will generate the final image with both left-eye and right-eye ones jointed together,
# with left-eye image on top and right-eye image on bottom.
# Recommended resolution for YouTube virtual reality videos is [3840, 2160]
cam.resolution = [500, 500]
cam.position = ds.arr([0.4, 0.5, 0.5], "code_length")
cam.switch_orientation(normal_vector=normal_vector, north_vector=north_vector)
# In (stereo)spherical camera, camera width is not used since the entire volume
# will be rendered
# Set the distance between left-eye and right-eye.
cam.lens.disparity = ds.domain_width[0] * 1.0e-3
sc.add_source(vol)
sc.save("lens_stereo-spherical.png", sigma_clip=6.0)

import yt

ds = yt.load('Enzo_64/DD0043/data0043')

sc = yt.create_scene(ds, lens_type='perspective')

source = sc[0]

source.set_field('density')
source.set_log(True)

# Set up the camera parameters: focus, width, resolution, and image orientation
sc.camera.focus = ds.domain_center
sc.camera.resolution = 1024
sc.camera.north_vector = [0, 0, 1]
sc.camera.position = [1.7, 1.7, 1.7]

# You may need to adjust the alpha values to get an image with good contrast.
# For the annotate_domain call, the fourth value in the color tuple is the
# alpha value.
sc.annotate_axes(alpha=.02)
sc.annotate_domain(ds, color=[1, 1, 1, .01])

text_string = "T = {} Gyr".format(float(ds.current_time.to('Gyr')))

# save an annotated version of the volume rendering including a representation
# of the transfer function and a nice label showing the simulation time.
sc.save_annotated("vol_annotated.png", sigma_clip=6,
                  text_annotate=[[(.1, 1.05), text_string]])

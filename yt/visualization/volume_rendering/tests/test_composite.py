from yt.mods import *
from yt.testing import \
    fake_random_pf
from yt.visualization.volume_rendering.scene import Scene
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.zbuffer_array import ZBuffer
from yt.visualization.volume_rendering.render_source import VolumeSource,\
    OpaqueSource
from yt.utilities.lib.misc_utilities import \
    lines
np.random.seed(0)

pf = fake_random_pf(64)
ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 3.)
pf.field_info[pf.field_list[0]].take_log=False

sc = Scene()
cam = Camera(ds)
sc.set_default_camera(cam)
vr = VolumeSource(ds, field=pf.field_list[0])
vr.transfer_function.clear()
vr.transfer_function.grey_opacity=True
vr.transfer_function.map_to_colormap(0.0, 1.0, scale=10.0, colormap="Reds")
sc.add_source(vr)

op = OpaqueSource()
empty = 0.0 * sc.default_camera.lens.new_image(cam)
z = np.ones(empty.shape[:2]) * np.inf

zbuff = ZBuffer(empty, z)
op.set_zbuffer(zbuff)

sc.add_source(op)

cam.set_width( pf.domain_width )

# DRAW SOME LINES
npoints = 100
vertices = 0.5 * np.random.random([npoints, 3])
#vertices = np.array([pf.domain_center, [2., 2., 1.]])
cam.lens.setup_box_properties(cam)
px, py, dz = cam.lens.project_to_plane(cam, vertices)
print dz
colors = np.random.random([npoints, 4])
colors[:,3] = 1.0 
lines(empty, px, py, colors, 24)
#empty[px, py, :] = 1.0
#z[px, py] = dz
dummy = -np.ones_like(empty)
lines(dummy, px, py, np.vstack([dz, dz]), 24)
print dummy[dummy!=-1]
z[:,:] = dummy[:,:,3]
z[z==-1] = np.inf
print z.min(), z.max()

im = sc.composite()
im = ImageArray(im.d)
im.write_png("composite.png")
#write_bitmap(zbuff.rgba[:, :, :3], 'composite.png')

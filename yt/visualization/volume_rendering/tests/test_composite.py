from yt.mods import *
from yt.testing import \
    fake_random_pf
from yt.visualization.volume_rendering.scene import Scene, RenderScene, \
    create_volume_rendering
from yt.visualization.volume_rendering.camera import Camera
from yt.visualization.volume_rendering.zbuffer_array import ZBuffer
from yt.visualization.volume_rendering.render_source import VolumeSource,\
    OpaqueSource


pf = fake_random_pf(64)
ds = pf.h.sphere(pf.domain_center, pf.domain_width[0] / 2)

sc = Scene()
cam = Camera(ds)
vr = VolumeSource(ds, field=('gas', 'density'))
sc.add_source(vr)
vr.build_defaults()

op = OpaqueSource()
op.set_scene(sc)
empty = 0.0 * vr.new_image()
z = np.ones(empty.shape[:2]) * np.inf
zbuff = ZBuffer(empty, z)
op.set_zbuffer(zbuff)

sc.add_source(op)

sc.render('test.png')

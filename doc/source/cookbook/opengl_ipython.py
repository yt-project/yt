
import yt
import numpy as np
from yt.visualization.volume_rendering.interactive_loop import \
    SceneGraph, BlockCollection, TrackballCamera, RenderingContext
import yt.visualization.volume_rendering.glfw_inputhook

rc = RenderingContext(1280, 960)

scene = SceneGraph()
collection = BlockCollection()

N = 64
x, y, z = np.mgrid[0:1:N*1j, 0:1:N*1j, 0:1:N*1j]
c = [-0.05, -0.05, -0.05]
oor2 = ((x-c[0])**2 + (y-c[1])**2 + (z-c[2])**2)**-0.5
np.clip(oor2, 1e-6, 1e60, oor2)
data = {'x_field': 10**x, 'y_field': 10**y, 'z_field': 10**z, 'sphere':oor2}

ds = yt.load_uniform_grid(data, [N, N, N], 1.0, nprocs=4)
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

dd = ds.all_data()
#collection.add_data(dd, "sphere")
collection.add_data(dd, "density")

scene.add_collection(collection)

position = (1.0, 1.0, 1.0)
c = TrackballCamera(position = position, focus = ds.domain_center,
                    near_plane = 0.1)

callbacks = rc.setup_loop(scene, c)
rl = rc(scene, c, callbacks)

# To make this work from IPython, you now do:
#
# yt.visualization.volume_rendering.glfw_inputhook.inputhook_manager.enable_gui("glfw", app=rl)

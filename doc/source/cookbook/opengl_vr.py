import yt
import numpy as np
from yt.visualization.volume_rendering.interactive_loop import \
    SceneGraph, BlockCollection, Camera, RenderingContext

rc = RenderingContext(1280, 960)

scene = SceneGraph()
collection = BlockCollection()

N = 64
x, y, z = np.mgrid[0:1:N*1j, 0:1:N*1j, 0:1:N*1j]
c = [-0.05, -0.05, -0.05]
oor2 = ((x-c[0])**2 + (y-c[1])**2 + (z-c[2])**2)**-0.5
np.clip(oor2, 1e-6, 1e60, oor2)
data = {'x_field': 10**x, 'y_field': 10**y, 'z_field': 10**z, 'sphere':oor2}

ds = yt.load_uniform_grid(data, [N, N, N], 1.0, nprocs=1)

dd = ds.all_data()
collection.add_data(dd, "z_field")

scene.add_collection(collection)

position = (-5.0, -5.0, 5.0)
c = Camera(position = position, focus = ds.domain_center)

rc.start_loop(scene, c)

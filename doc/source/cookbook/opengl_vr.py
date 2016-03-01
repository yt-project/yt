import yt
import numpy as np
from yt.visualization.volume_rendering.interactive_loop import \
    SceneGraph, BlockCollection, TrackballCamera, RenderingContext

rc = RenderingContext(1280, 960)

scene = SceneGraph()
collection = BlockCollection()

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

dd = ds.all_data()
collection.add_data(dd, "density")

scene.add_collection(collection)

position = (1.0, 1.0, 1.0)
c = TrackballCamera(position = position, focus = ds.domain_center,
                    near_plane = 0.1)

rc.start_loop(scene, c)

import yt
from yt.visualization.volume_rendering.interactive_vr import \
    SceneGraph, BlockCollection, TrackballCamera
from yt.visualization.volume_rendering.interactive_loop import \
    RenderingContext
from yt.visualization.volume_rendering import glfw_inputhook 

rc = RenderingContext(1280, 960)

scene = SceneGraph()
collection = BlockCollection()

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

dd = ds.all_data()
collection.add_data(dd, "density")

scene.add_collection(collection)

position = (1.0, 1.0, 1.0)
c = TrackballCamera(position=position, focus=ds.domain_center,
                    near_plane=0.1)

callbacks = rc.setup_loop(scene, c)
rl = rc(scene, c, callbacks)

# To make this work from IPython execute:
#
# glfw_inputhook.inputhook_manager.enable_gui("glfw", app=rl)

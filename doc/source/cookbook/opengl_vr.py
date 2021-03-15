import yt
from yt.visualization.volume_rendering.interactive_loop import RenderingContext
from yt.visualization.volume_rendering.interactive_vr import (
    BlockCollection,
    SceneGraph,
    TrackballCamera,
)

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create GLUT window
rc = RenderingContext(1280, 960)

# Create a 3d Texture from all_data()
collection = BlockCollection()
dd = ds.all_data()
collection.add_data(dd, "density")

# Initialize basic Scene and pass the data
scene = SceneGraph()
scene.add_collection(collection)

# Create default camera
position = (1.0, 1.0, 1.0)
c = TrackballCamera(position=position, focus=ds.domain_center, near_plane=0.1)

# Start rendering loop
rc.start_loop(scene, c)

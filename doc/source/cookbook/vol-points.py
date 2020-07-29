import numpy as np

import yt
from yt.units import kpc
from yt.visualization.volume_rendering.api import PointSource

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

sc = yt.create_scene(ds)

np.random.seed(1234567)

npoints = 1000

# Random particle positions
vertices = np.random.random([npoints, 3]) * 200 * kpc

# Random colors
colors = np.random.random([npoints, 4])

# Set alpha value to something that produces a good contrast with the volume
# rendering
colors[:, 3] = 0.1

points = PointSource(vertices, colors=colors)
sc.add_source(points)

sc.camera.width = 300 * kpc

sc.save(sigma_clip=5)

import numpy as np

import yt
from yt.units import kpc
from yt.visualization.volume_rendering.api import LineSource

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

sc = yt.create_scene(ds)

np.random.seed(1234567)

nlines = 50
vertices = (np.random.random([nlines, 2, 3]) - 0.5) * 200 * kpc
colors = np.random.random([nlines, 4])
colors[:, 3] = 0.1

lines = LineSource(vertices, colors)
sc.add_source(lines)

sc.camera.width = 300 * kpc

sc.save(sigma_clip=4.0)

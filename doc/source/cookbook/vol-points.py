import yt
import numpy as np
from yt.visualization.volume_rendering.api import PointSource
from yt.units import kpc

ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')

sc = yt.create_scene(ds)

np.random.seed(1234567)

npoints = 1000
vertices = np.random.random([npoints, 3])*200*kpc
colors = np.random.random([npoints, 4])
colors[:, 3] = 1
colors *= 0.1

points = PointSource(vertices, colors=colors)
sc.add_source(points)

sc.camera.width = 300*kpc

sc.save()

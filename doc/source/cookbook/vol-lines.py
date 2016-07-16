import yt
import numpy as np
from yt.visualization.volume_rendering.api import LineSource
from yt.units import kpc

np.random.seed(1234567)

ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')

im, sc = yt.volume_render(ds)

nlines = 50
vertices = (np.random.random([nlines, 2, 3]) - 0.5) * 250 * kpc
colors = np.random.random([nlines, 4])
colors[:, 3] = 1.0

lines = LineSource(vertices, colors)
sc.add_source(lines)

sc.camera.width = 300*kpc

sc.save()

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import yt

# Load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere object centered on the highest density point in the simulation
# with radius 1 Mpc
sphere = ds.sphere("max", (1.0, "Mpc"))

# Identify the isodensity surface in this sphere with density = 1e-24 g/cm^3
surface = ds.surface(sphere, ("gas", "density"), 1e-24)

# Color this isodensity surface according to the log of the temperature field
colors = yt.apply_colormap(np.log10(surface[("gas", "temperature")]), cmap_name="hot")

# Create a 3D matplotlib figure for visualizing the surface
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

# Set the surface colors in the right scaling [0,1]
p3dc.set_facecolors(colors[0, :, :] / 255.0)
ax.add_collection(p3dc)

# Let's keep the axis ratio fixed in all directions by taking the maximum
# extent in one dimension and make it the bounds in all dimensions
max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
bounds = np.zeros([3, 2])
bounds[:, 0] = centers[:] - max_extent / 2
bounds[:, 1] = centers[:] + max_extent / 2
ax.auto_scale_xyz(bounds[0, :], bounds[1, :], bounds[2, :])

# Save the figure
plt.savefig(f"{ds}_Surface.png")

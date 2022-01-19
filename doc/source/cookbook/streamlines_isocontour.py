import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import yt
from yt.visualization.api import Streamlines

# Load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Define c: the center of the box, N: the number of streamlines,
# scale: the spatial scale of the streamlines relative to the boxsize,
# and then pos: the random positions of the streamlines.
c = ds.arr([0.5] * 3, "code_length")
N = 30
scale = ds.quan(15, "kpc").in_units("code_length")  # 15 kpc in code units
pos_dx = np.random.random((N, 3)) * scale - scale / 2.0
pos = c + pos_dx

# Create the streamlines from these positions with the velocity fields as the
# fields to be traced
streamlines = Streamlines(
    ds,
    pos,
    ("gas", "velocity_x"),
    ("gas", "velocity_y"),
    ("gas", "velocity_z"),
    length=1.0,
)
streamlines.integrate_through_volume()

# Create a 3D matplotlib figure for visualizing the streamlines
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

# Trace the streamlines through the volume of the 3D figure
for stream in streamlines.streamlines:
    stream = stream[np.all(stream != 0.0, axis=1)]

    # Make the colors of each stream vary continuously from blue to red
    # from low-x to high-x of the stream start position (each color is R, G, B)
    # can omit and just set streamline colors to a fixed color
    x_start_pos = ds.arr(stream[0, 0], "code_length")
    x_start_pos -= ds.arr(0.5, "code_length")
    x_start_pos /= scale
    x_start_pos += 0.5
    color = np.array([x_start_pos, 0, 1 - x_start_pos])

    # Plot the stream in 3D
    ax.plot3D(stream[:, 0], stream[:, 1], stream[:, 2], alpha=0.3, color=color)

# Create a sphere object centered on the highest density point in the simulation
# with radius = 1 Mpc
sphere = ds.sphere("max", (1.0, "Mpc"))

# Identify the isodensity surface in this sphere with density = 1e-24 g/cm^3
surface = ds.surface(sphere, ("gas", "density"), 1e-24)

# Color this isodensity surface according to the log of the temperature field
colors = yt.apply_colormap(np.log10(surface[("gas", "temperature")]), cmap_name="hot")

# Render this surface
p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
colors = colors[0, :, :] / 255.0  # scale to [0,1]
colors[:, 3] = 0.3  # alpha = 0.3
p3dc.set_facecolors(colors)
ax.add_collection(p3dc)

# Save the figure
plt.savefig("streamlines_isocontour.png")

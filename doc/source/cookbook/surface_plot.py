### THIS RECIPE IS CURRENTLY BROKEN IN YT-3.0
### DO NOT TRUST THIS RECIPE UNTIL THIS LINE IS REMOVED

import yt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

# Load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere object centered on the highest density point in the simulation
# with radius 1 Mpc
sphere = ds.sphere("max", (1.0, "Mpc"))

# Identify the isodensity surface in this sphere with density = 1e-24 g/cm^3
surface = ds.surface(sphere, "density", 1e-24)

# Color this isodensity surface according to the log of the temperature field
colors = yt.apply_colormap(np.log10(surface["temperature"]), cmap_name="hot")

# Create a 3D matplotlib figure for visualizing the surface
fig = plt.figure()
ax = fig.gca(projection='3d')
p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

# Set the surface colors in the right scaling [0,1]
p3dc.set_facecolors(colors[0,:,:]/255.)
ax.add_collection(p3dc)
ax.auto_scale_xyz(surface.vertices[0,:], surface.vertices[1,:], surface.vertices[2,:])
ax.set_aspect(1.0)

# Save the figure
plt.savefig("%s_Surface.png" % ds)

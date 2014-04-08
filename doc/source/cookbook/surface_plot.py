from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from yt.mods import *

pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
sphere = pf.sphere("max", (1.0, "mpc"))
surface = pf.surface(sphere, "density", 1e-25)
colors = apply_colormap(np.log10(surface["temperature"]), cmap_name="hot")

fig = plt.figure()
ax = fig.gca(projection='3d')
p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

p3dc.set_facecolors(colors[0,:,:]/255.)
ax.add_collection(p3dc)
ax.auto_scale_xyz(surface.vertices[0,:], surface.vertices[1,:], surface.vertices[2,:])
ax.set_aspect(1.0)

plt.savefig("%s_Surface.png" % pf)

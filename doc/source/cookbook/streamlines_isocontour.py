from yt.mods import *
from yt.visualization.api import Streamlines

pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
c = np.array([0.5]*3)
N = 30
scale = 15.0/pf['kpc']
pos_dx = np.random.random((N,3))*scale-scale/2.
pos = c+pos_dx

streamlines = Streamlines(pf,pos,'velocity_x', 'velocity_y', 'velocity_z', length=1.0) 
streamlines.integrate_through_volume()

import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
fig=pl.figure() 
ax = Axes3D(fig)
for stream in streamlines.streamlines:
    stream = stream[np.all(stream != 0.0, axis=1)]
    ax.plot3D(stream[:,0], stream[:,1], stream[:,2], alpha=0.1)


sphere = pf.sphere("max", (1.0, "mpc"))
surface = pf.surface(sphere, "density", 1e-24)
colors = apply_colormap(np.log10(surface["temperature"]), cmap_name="hot")

p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
colors = colors[0,:,:]/255.
colors[:,3] = 0.3
p3dc.set_facecolors(colors)
ax.add_collection(p3dc)

pl.savefig('streamlines_isocontour.png')


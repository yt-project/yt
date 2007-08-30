#
# How to use both yt and Matplotlib to make a 3x3 plot collection,
# along all three axes and with three different widths.
# (For each width, the three axes have synced colorbars.)
#

fn = "DataDump0043.dir/DataDump0043"
field = "Temperature"
cmap = "hot"
widths = [(1,"pc"),(1000,"au"),(10,"au")]

# Below this line, you shouldn't need to change anything if you just want the
# vanilla plot

import yt.lagos as lagos
import yt.raven as raven
import matplotlib.figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
# backend_agg is offscreen

a = lagos.EnzoStaticOutput(fn)
v, c = a.hierarchy.findMax("Density") # Center on max
# c == center

fig = matplotlib.figure.Figure((9,9))
canvas = FigureCanvasAgg(fig)
fig.subplots_adjust(hspace=0,wspace=0,bottom=0.0, top=1.0, left=0.0, right=1.0)

ax = [fig.add_subplot(3,3,i+1) for i in range(9)] # Set up our nine 'axes'

slices = []
slices.append(a.h.slice(axis=0, fields=field, coord=c[0], center=c))
slices.append(a.h.slice(axis=1, fields=field, coord=c[1], center=c))
slices.append(a.h.slice(axis=2, fields=field, coord=c[2], center=c))

# This could also be done with list comprehensions thusly:
#
# slices = [a.h.slice(ax, "Temperature", coord=c[ax], center=c)
#               for ax in [0,1,2] ]

# We set up a PlotCollection to contain each axis
pc = []

for i in range(3):
    pc.append(raven.PlotCollection(a))
    for j in range(3):
        pc[-1].addPlot(raven.be.SlicePlot(proj[j], field, fig, ax[i*3+j], False))
        pc[-1].plots[-1].set_cmap(cmap)
    w, u = widths[i]
    pc[-1].set_width(w,u) # This propagates to all the plots

for p in pc:
    vmin = 1e30
    vmax = 1e-30
    for j in range(3):
        vmin = min(vmin, p.plots[-1].norm.vmin)
        vmax = max(vmax, p.plots[-1].norm.vmax)
    p.set_zlim(vmin,vmax)

# Note that we're calling savefig on the figure, which makes it 
# a matplotlib call.

fig.savefig("ThreeByThreePlot.png")

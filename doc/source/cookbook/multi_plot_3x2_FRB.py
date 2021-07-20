import numpy as np
from matplotlib.colors import LogNorm

import yt
from yt.visualization.api import get_multi_plot

fn = "Enzo_64/RD0006/RedshiftOutput0006"  # dataset to load

# load data and get center value and center location as maximum density location
ds = yt.load(fn)
v, c = ds.find_max(("gas", "density"))

# set up our Fixed Resolution Buffer parameters: a width, resolution, and center
width = (1.0, "unitary")
res = [1000, 1000]
#  get_multi_plot returns a containing figure, a list-of-lists of axes
#   into which we can place plots, and some axes that we'll put
#   colorbars.

#  it accepts: # of x-axis plots, # of y-axis plots, and how the
#  colorbars are oriented (this also determines where they go: below
#  in the case of 'horizontal', on the right in the case of
#  'vertical'), bw is the base-width in inches (4 is about right for
#  most cases)

orient = "horizontal"
fig, axes, colorbars = get_multi_plot(2, 3, colorbar=orient, bw=6)

# Now we follow the method of "multi_plot.py" but we're going to iterate
# over the columns, which will become axes of slicing.
plots = []
for ax in range(3):
    sli = ds.slice(ax, c[ax])
    frb = sli.to_frb(width, res)
    den_axis = axes[ax][0]
    temp_axis = axes[ax][1]

    # here, we turn off the axes labels and ticks, but you could
    # customize further.
    for ax in (den_axis, temp_axis):
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

    # converting our fixed resolution buffers to NDarray so matplotlib can
    # render them
    dens = np.array(frb[("gas", "density")])
    temp = np.array(frb[("gas", "temperature")])

    plots.append(den_axis.imshow(dens, norm=LogNorm()))
    plots[-1].set_clim((5e-32, 1e-29))
    plots[-1].set_cmap("bds_highcontrast")

    plots.append(temp_axis.imshow(temp, norm=LogNorm()))
    plots[-1].set_clim((1e3, 1e8))
    plots[-1].set_cmap("hot")

# Each 'cax' is a colorbar-container, into which we'll put a colorbar.
# the zip command creates triples from each element of the three lists
# .  Note that it cuts off after the shortest iterator is exhausted,
# in this case, titles.
titles = [
    r"$\mathrm{density}\ (\mathrm{g\ cm^{-3}})$",
    r"$\mathrm{temperature}\ (\mathrm{K})$",
]
for p, cax, t in zip(plots, colorbars, titles):
    # Now we make a colorbar, using the 'image' we stored in plots
    # above. note this is what is *returned* by the imshow method of
    # the plots.
    cbar = fig.colorbar(p, cax=cax, orientation=orient)
    cbar.set_label(t)

# And now we're done!
fig.savefig(f"{ds}_3x2.png")

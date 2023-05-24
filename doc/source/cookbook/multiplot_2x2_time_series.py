import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

import yt

fns = [
    "Enzo_64/DD0005/data0005",
    "Enzo_64/DD0015/data0015",
    "Enzo_64/DD0025/data0025",
    "Enzo_64/DD0035/data0035",
]

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot with a single
# shared narrow colorbar on the right hand side of the multipanel plot. Axes
# labels are drawn for all plots since we're slicing along different directions
# for each plot.
grid = AxesGrid(
    fig,
    (0.075, 0.075, 0.85, 0.85),
    nrows_ncols=(2, 2),
    axes_pad=0.05,
    label_mode="L",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)

for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn)  # load data
    p = yt.ProjectionPlot(ds, "z", ("gas", "density"), width=(55, "Mpccm"))

    # Ensure the colorbar limits match for all plots
    p.set_zlim(("gas", "density"), 1e-4, 1e-2)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", "density")]
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p.render()

plt.savefig("multiplot_2x2_time_series.png")

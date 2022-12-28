import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

import yt

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
grid = AxesGrid(
    fig,
    (0.085, 0.085, 0.83, 0.83),
    nrows_ncols=(1, 2),
    axes_pad=0.05,
    label_mode="L",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
    aspect=False,
)

for i, SnapNum in enumerate([10, 40]):
    # Load the data and create a single plot
    ds = yt.load("enzo_tiny_cosmology/DD00%2d/DD00%2d" % (SnapNum, SnapNum))
    ad = ds.all_data()
    p = yt.PhasePlot(
        ad,
        ("gas", "density"),
        ("gas", "temperature"),
        [
            ("gas", "mass"),
        ],
        weight_field=None,
    )

    # Ensure the axes and colorbar limits match for all plots
    p.set_xlim(1.0e-32, 8.0e-26)
    p.set_ylim(1.0e1, 2.0e7)
    p.set_zlim(("gas", "mass"), 1e42, 1e46)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", "mass")]
    plot.figure = fig
    plot.axes = grid[i].axes
    if i == 0:
        plot.cax = grid.cbar_axes[i]

    # Actually redraws the plot.
    p.render()

    # Modify the axes properties **after** p.render() so that they
    # are not overwritten.
    plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))

plt.savefig("multiplot_phaseplot.png")

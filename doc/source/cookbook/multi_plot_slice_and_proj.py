import numpy as np
from matplotlib.colors import LogNorm

import yt
from yt.visualization.base_plot_types import get_multi_plot

fn = "GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150"  # dataset to load
orient = "horizontal"

ds = yt.load(fn)  # load data

# There's a lot in here:
#   From this we get a containing figure, a list-of-lists of axes into which we
#   can place plots, and some axes that we'll put colorbars.
# We feed it:
#   Number of plots on the x-axis, number of plots on the y-axis, and how we
#   want our colorbars oriented.  (This governs where they will go, too.
#   bw is the base-width in inches, but 4 is about right for most cases.
fig, axes, colorbars = get_multi_plot(3, 2, colorbar=orient, bw=4)

slc = yt.SlicePlot(
    ds,
    "z",
    fields=[("gas", "density"), ("gas", "temperature"), ("gas", "velocity_magnitude")],
)
proj = yt.ProjectionPlot(ds, "z", ("gas", "density"), weight_field=("gas", "density"))

slc_frb = slc.data_source.to_frb((1.0, "Mpc"), 512)
proj_frb = proj.data_source.to_frb((1.0, "Mpc"), 512)

dens_axes = [axes[0][0], axes[1][0]]
temp_axes = [axes[0][1], axes[1][1]]
vels_axes = [axes[0][2], axes[1][2]]

for dax, tax, vax in zip(dens_axes, temp_axes, vels_axes):

    dax.xaxis.set_visible(False)
    dax.yaxis.set_visible(False)
    tax.xaxis.set_visible(False)
    tax.yaxis.set_visible(False)
    vax.xaxis.set_visible(False)
    vax.yaxis.set_visible(False)

# Converting our Fixed Resolution Buffers to numpy arrays so that matplotlib
# can render them

slc_dens = np.array(slc_frb[("gas", "density")])
proj_dens = np.array(proj_frb[("gas", "density")])
slc_temp = np.array(slc_frb[("gas", "temperature")])
proj_temp = np.array(proj_frb[("gas", "temperature")])
slc_vel = np.array(slc_frb[("gas", "velocity_magnitude")])
proj_vel = np.array(proj_frb[("gas", "velocity_magnitude")])

plots = [
    dens_axes[0].imshow(slc_dens, origin="lower", norm=LogNorm()),
    dens_axes[1].imshow(proj_dens, origin="lower", norm=LogNorm()),
    temp_axes[0].imshow(slc_temp, origin="lower"),
    temp_axes[1].imshow(proj_temp, origin="lower"),
    vels_axes[0].imshow(slc_vel, origin="lower", norm=LogNorm()),
    vels_axes[1].imshow(proj_vel, origin="lower", norm=LogNorm()),
]

plots[0].set_clim((1.0e-27, 1.0e-25))
plots[0].set_cmap("bds_highcontrast")
plots[1].set_clim((1.0e-27, 1.0e-25))
plots[1].set_cmap("bds_highcontrast")
plots[2].set_clim((1.0e7, 1.0e8))
plots[2].set_cmap("hot")
plots[3].set_clim((1.0e7, 1.0e8))
plots[3].set_cmap("hot")
plots[4].set_clim((1e6, 1e8))
plots[4].set_cmap("gist_rainbow")
plots[5].set_clim((1e6, 1e8))
plots[5].set_cmap("gist_rainbow")

titles = [
    r"$\mathrm{Density}\ (\mathrm{g\ cm^{-3}})$",
    r"$\mathrm{Temperature}\ (\mathrm{K})$",
    r"$\mathrm{Velocity Magnitude}\ (\mathrm{cm\ s^{-1}})$",
]

for p, cax, t in zip(plots[0:6:2], colorbars, titles):
    cbar = fig.colorbar(p, cax=cax, orientation=orient)
    cbar.set_label(t)

# And now we're done!
fig.savefig(f"{ds}_3x2")

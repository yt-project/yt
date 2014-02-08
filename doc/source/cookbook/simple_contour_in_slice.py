from yt.mods import *

# Load the data file.
pf = load("Sedov_3d/sedov_hdf5_chk_0002")

# Make a traditional slice plot.
sp = SlicePlot(pf,"x","dens")

# Overlay the slice plot with thick red contours of density.
sp.annotate_contour("dens", ncont=3, clim=(1e-2,1e-1), label=True,
                    plot_args={"colors": "red",
                               "linewidths": 2})

# What about some nice temperature contours in blue?
sp.annotate_contour("temp", ncont=3, clim=(1e-8,1e-6), label=True,
                    plot_args={"colors": "blue",
                               "linewidths": 2})

# This is the plot object.
po = sp.plots["dens"]

# Turn off the colormap image, leaving just the contours.
po.axes.images[0].set_visible(False)

# Remove the colorbar and its label.
po.figure.delaxes(po.figure.axes[1])

# Save it and ask for a close fit to get rid of the space used by the colorbar.
sp.save(mpl_kwargs={'bbox_inches':'tight'})

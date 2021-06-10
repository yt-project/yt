import numpy as np

import yt

# Load the dataset.
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create a slice object
slc = yt.SlicePlot(ds, "x", ("gas", "density"), width=(800.0, "kpc"))

# Get a reference to the matplotlib axes object for the plot
ax = slc.plots[("gas", "density")].axes

# Let's adjust the x axis tick labels
for label in ax.xaxis.get_ticklabels():
    label.set_color("red")
    label.set_fontsize(16)

# Get a reference to the matplotlib figure object for the plot
fig = slc.plots[("gas", "density")].figure

# And create a mini-panel of a gaussian histogram inside the plot
rect = (0.2, 0.2, 0.2, 0.2)
new_ax = fig.add_axes(rect)

n, bins, patches = new_ax.hist(
    np.random.randn(1000) + 20, 50, facecolor="black", edgecolor="black"
)

# Make sure its visible
new_ax.tick_params(colors="white")

# And label it
la = new_ax.set_xlabel("Dinosaurs per furlong")
la.set_color("white")

slc.save()

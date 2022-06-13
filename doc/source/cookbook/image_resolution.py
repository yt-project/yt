import numpy as np

import yt

# Load the dataset.  We'll work with a some Gadget data to illustrate all
# the different ways in which the effective resolution can vary.  Specifically,
# we'll use the GadgetDiskGalaxy dataset available at
#  http://yt-project.org/data/GadgetDiskGalaxy.tar.gz

# load the data with a refinement criteria of 2 particle per cell
# n.b. -- in yt-4.0, n_ref no longer exists as the data is no longer
# deposited only a grid.  At present (03/15/2019), there is no way to
# handle non-gas data in Gadget snapshots, though that is work in progress
if int(yt.__version__[0]) < 4:
    # increasing n_ref will result in a "lower resolution" (but faster) image,
    # while decreasing it will go the opposite way
    ds = yt.load("GadgetDiskGalaxy/snapshot_200.hdf5", n_ref=16)
else:
    ds = yt.load("GadgetDiskGalaxy/snapshot_200.hdf5")

# Create projections of the density (max value in each resolution element in the image):
prj = yt.ProjectionPlot(
    ds, "x", ("gas", "density"), method="max", center="max", width=(100, "kpc")
)

# nicen up the plot by using a better interpolation:
plot = prj.plots[list(prj.plots)[0]]
ax = plot.axes
img = ax.images[0]
img.set_interpolation("bicubic")

# nicen up the plot by setting the background color to the minimum of the colorbar
prj.set_background_color(("gas", "density"))

# vary the buff_size -- the number of resolution elements in the actual visualization
# set it to 2000x2000
buff_size = 2000
prj.set_buff_size(buff_size)

# set the figure size in inches
figure_size = 10
prj.set_figure_size(figure_size)

# if the image does not fill the plot (as is default, since the axes and
# colorbar contribute as well), then figuring out the proper dpi for a given
# buff_size and figure_size is non-trivial -- it requires finding the bbox
# for the actual image:
bounding_box = ax.get_position()
# we're going to scale to the larger of the two sides
image_size = figure_size * max([bounding_box.width, bounding_box.height])
# now save with a dpi that's scaled to the buff_size:
dpi = np.rint(np.ceil(buff_size / image_size))
prj.save("with_axes_colorbar.png", mpl_kwargs=dict(dpi=dpi))

# in the case where the image fills the entire plot (i.e. if the axes and colorbar
# are turned off), it's trivial to figure out the correct dpi from the buff_size and
# figure_size (or vice versa):

# hide the colorbar:
prj.hide_colorbar()

# hide the axes, while still keeping the background color correct:
prj.hide_axes(draw_frame=True)

# save with a dpi that makes sense:
dpi = np.rint(np.ceil(buff_size / figure_size))
prj.save("no_axes_colorbar.png", mpl_kwargs=dict(dpi=dpi))

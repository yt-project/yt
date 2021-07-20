import yt

fn = "IsolatedGalaxy/galaxy0030/galaxy0030"  # dataset to load

ds = yt.load(fn)  # load data

# Now we need a center of our volume to render.  Here we'll just use
# 0.5,0.5,0.5, because volume renderings are not periodic.
c = [0.5, 0.5, 0.5]

# Our image plane will be normal to some vector.  For things like collapsing
# objects, you could set it the way you would a cutting plane -- but for this
# dataset, we'll just choose an off-axis value at random.  This gets normalized
# automatically.
L = [0.5, 0.4, 0.7]

# Our "width" is the width of the image plane as well as the depth.
# The first element is the left to right width, the second is the
# top-bottom width, and the last element is the back-to-front width
# (all in code units)
W = [0.04, 0.04, 0.4]

# The number of pixels along one side of the image.
# The final image will have Npixel^2 pixels.
Npixels = 512

# Now we call the off_axis_projection function, which handles the rest.
# Note that we set no_ghost equal to False, so that we *do* include ghost
# zones in our data.  This takes longer to calculate, but the results look
# much cleaner than when you ignore the ghost zones.
# Also note that we set the field which we want to project as "density", but
# really we could use any arbitrary field like "temperature", "metallicity"
# or whatever.
image = yt.off_axis_projection(ds, c, L, W, Npixels, ("gas", "density"), no_ghost=False)

# Image is now an NxN array representing the intensities of the various pixels.
# And now, we call our direct image saver.  We save the log of the result.
yt.write_projection(
    image,
    "offaxis_projection_colorbar.png",
    colorbar_label="Column Density (cm$^{-2}$)",
)

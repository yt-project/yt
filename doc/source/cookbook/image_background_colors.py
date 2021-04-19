import yt

# This shows how to save ImageArray objects, such as those returned from
# volume renderings, to pngs with varying backgrounds.

# First we use the simple_volume_rendering.py recipe from above to generate
# a standard volume rendering.

ds = yt.load("Enzo_64/DD0043/data0043")
im, sc = yt.volume_render(ds, ("gas", "density"))
im.write_png("original.png", sigma_clip=8.0)

# Our image array can now be transformed to include different background
# colors.  By default, the background color is black.  The following
# modifications can be used on any image array.

# write_png accepts a background keyword argument that defaults to 'black'.
# Other choices include:
# black (0.,0.,0.,1.)
# white (1.,1.,1.,1.)
# None  (0.,0.,0.,0.) <-- Transparent!
# any rgba list/array: [r,g,b,a], bounded by 0..1

# We include the sigma_clip=8 keyword here to bring out more contrast between
# the background and foreground, but it is entirely optional.

im.write_png("black_bg.png", background="black", sigma_clip=8.0)
im.write_png("white_bg.png", background="white", sigma_clip=8.0)
im.write_png("green_bg.png", background=[0.0, 1.0, 0.0, 1.0], sigma_clip=8.0)
im.write_png("transparent_bg.png", background=None, sigma_clip=8.0)

import yt

# Load the dataset.
ds = yt.load("enzo_tiny_cosmology/RD0009/RD0009")

# Create a volume rendering, which will determine data bounds, use the first
# acceptable field in the field_list, and set up a default transfer function.

# Render and save output images with different levels of sigma clipping.
# Sigma clipping removes the highest intensity pixels in a volume render,
# which affects the overall contrast of the image.
sc = yt.create_scene(ds, field=("gas", "density"))
sc.save("clip_0.png")
sc.save("clip_2.png", sigma_clip=2)
sc.save("clip_4.png", sigma_clip=4)
sc.save("clip_6.png", sigma_clip=6)

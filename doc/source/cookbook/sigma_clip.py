import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a volume rendering, which will determine data bounds, use the first
# acceptable field in the field_list, and set up a default transfer function.

# Render and save output images with different levels of sigma clipping
sc = yt.create_scene(ds, field=('gas', 'density'))
sc.render()
sc.save('raw.png')
sc.save('clip_2.png', sigma_clip=2)
sc.save('clip_4.png', sigma_clip=4)
sc.save('clip_6.png', sigma_clip=6)

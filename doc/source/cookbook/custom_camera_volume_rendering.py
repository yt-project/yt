import yt

# Load the dataset
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a volume rendering
sc = yt.create_scene(ds, field=("gas", "density"))

# Now increase the resolution
sc.camera.resolution = (1024, 1024)

# Set the camera focus to a position that is offset from the center of
# the domain
sc.camera.focus = ds.arr([0.3, 0.3, 0.3], "unitary")

# Move the camera position to the other side of the dataset
sc.camera.position = ds.arr([0, 0, 0], "unitary")

# save to disk with a custom filename and apply sigma clipping to eliminate
# very bright pixels, producing an image with better contrast.
sc.render()
sc.save("custom.png", sigma_clip=4)

import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a 1D profile within a sphere of radius 100 kpc
# of the average temperature and average velocity_x
# vs. density, weighted by mass.
sphere = ds.sphere("c", (100.0, "kpc"))
plot = yt.ProfilePlot(
    sphere,
    ("gas", "density"),
    [("gas", "temperature"), ("gas", "velocity_x")],
    weight_field=("gas", "mass"),
)
plot.set_log(("gas", "velocity_x"), False)

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()

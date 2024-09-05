import yt

# load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# This is an object that describes the entire box
ad = ds.all_data()

# We plot the average velocity magnitude (mass-weighted) in our object
# as a function of density and temperature
plot = yt.PhasePlot(
    ad, ("gas", "density"), ("gas", "temperature"), ("gas", "velocity_magnitude")
)

# save the plot
plot.save()

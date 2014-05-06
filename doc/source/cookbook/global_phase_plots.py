import yt

# load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# This is an object that describes the entire box
ad = ds.h.all_data()

# We plot the average VelocityMagnitude (mass-weighted) in our object
# as a function of Density and temperature
plot = yt.PhasePlot(ad, "density", "temperature", "velocity_magnitude")

# save the plot
plot.save()

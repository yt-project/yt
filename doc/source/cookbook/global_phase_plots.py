from yt.mods import * # set up our namespace

# load the dataset
pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

# This is an object that describes the entire box
ad = pf.h.all_data()

# We plot the average VelocityMagnitude (mass-weighted) in our object 
# as a function of Density and Temperature
plot = PhasePlot(ad, "Density","Temperature","VelocityMagnitude")

# save the plot
plot.save()

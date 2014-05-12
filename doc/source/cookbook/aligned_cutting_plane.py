import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a 1 kpc radius sphere, centered on the maximum gas density.  Note
# that this sphere is very small compared to the size of our final plot,
# and it has a non-axially aligned L vector.
sp = ds.sphere("m", (1.0, "kpc"))

# Get the angular momentum vector for the sphere.
L = sp.quantities.angular_momentum_vector()

print "Angular momentum vector: {0}".format(L)

# Create an OffAxisSlicePlot on the object with the L vector as its normal
p = yt.OffAxisSlicePlot(ds, L, "density", sp.center, (15, "kpc"))
p.save()

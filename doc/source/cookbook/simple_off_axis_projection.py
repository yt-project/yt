from yt.mods import *

# Load the dataset.
pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a 1 kpc radius sphere, centered on the max density.  Note that this
# sphere is very small compared to the size of our final plot, and it has a
# non-axially aligned L vector.
sp = pf.sphere("center", (15.0, "kpc"))

# Get the angular momentum vector for the sphere.
L = sp.quantities["AngularMomentumVector"]()

print "Angular momentum vector: %s" % (L)

# Create an OffAxisSlicePlot on the object with the L vector as its normal
p = OffAxisProjectionPlot(pf, L, "density", sp.center, (25, "kpc"))
p.save()

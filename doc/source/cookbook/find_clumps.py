import numpy as np

import yt
from yt.data_objects.level_sets.api import Clump, find_clumps

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

data_source = ds.disk([0.5, 0.5, 0.5], [0.0, 0.0, 1.0], (8, "kpc"), (1, "kpc"))

# the field to be used for contouring
field = ("gas", "density")

# This is the multiplicative interval between contours.
step = 2.0

# Now we set some sane min/max values between which we want to find contours.
# This is how we tell the clump finder what to look for -- it won't look for
# contours connected below or above these threshold values.
c_min = 10 ** np.floor(np.log10(data_source[field]).min())
c_max = 10 ** np.floor(np.log10(data_source[field]).max() + 1)

# Now find get our 'base' clump -- this one just covers the whole domain.
master_clump = Clump(data_source, field)

# Add a "validator" to weed out clumps with less than 20 cells.
# As many validators can be added as you want.
master_clump.add_validator("min_cells", 20)

# Calculate center of mass for all clumps.
master_clump.add_info_item("center_of_mass")

# Begin clump finding.
find_clumps(master_clump, c_min, c_max, step)

# Save the clump tree as a reloadable dataset
fn = master_clump.save_as_dataset(fields=[("gas", "density"), ("all", "particle_mass")])

# We can traverse the clump hierarchy to get a list of all of the 'leaf' clumps
leaf_clumps = master_clump.leaves

# Get total cell and particle masses for each leaf clump
leaf_masses = [leaf.quantities.total_mass() for leaf in leaf_clumps]

# If you'd like to visualize these clumps, a list of clumps can be supplied to
# the "clumps" callback on a plot.  First, we create a projection plot:
prj = yt.ProjectionPlot(ds, 2, field, center="c", width=(20, "kpc"))

# Next we annotate the plot with contours on the borders of the clumps
prj.annotate_clumps(leaf_clumps)

# Save the plot to disk.
prj.save("clumps")

# Reload the clump dataset.
cds = yt.load(fn)

# Clump annotation can also be done with the reloaded clump dataset.

# Remove the original clump annotation
prj.clear_annotations()

# Get the leaves and add the callback.
leaf_clumps_reloaded = cds.leaves
prj.annotate_clumps(leaf_clumps_reloaded)
prj.save("clumps_reloaded")

# Query fields for clumps in the tree.
print(cds.tree["clump", "center_of_mass"])
print(cds.tree.children[0]["grid", "density"])
print(cds.tree.children[1]["all", "particle_mass"])

# Get all of the leaf clumps.
print(cds.leaves)
print(cds.leaves[0]["clump", "cell_mass"])

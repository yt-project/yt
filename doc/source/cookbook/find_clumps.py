import numpy as np

import yt
from yt.analysis_modules.level_sets.api import *

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

data_source = ds.disk([0.5, 0.5, 0.5], [0., 0., 1.],
                      (8, 'kpc'), (1, 'kpc'))

# the field to be used for contouring
field = ("gas", "density")

# This is the multiplicative interval between contours.
step = 2.0

# Now we set some sane min/max values between which we want to find contours.
# This is how we tell the clump finder what to look for -- it won't look for
# contours connected below or above these threshold values.
c_min = 10**np.floor(np.log10(data_source[field]).min()  )
c_max = 10**np.floor(np.log10(data_source[field]).max()+1)

# Now find get our 'base' clump -- this one just covers the whole domain.
master_clump = Clump(data_source, field)

# Add a "validator" to weed out clumps with less than 20 cells.
# As many validators can be added as you want.
master_clump.add_validator("min_cells", 20)

# Begin clump finding.
find_clumps(master_clump, c_min, c_max, step)

# Write out the full clump hierarchy.
write_clump_index(master_clump, 0, "%s_clump_hierarchy.txt" % ds)

# Write out only the leaf nodes of the hierarchy.
write_clumps(master_clump,0, "%s_clumps.txt" % ds)

# We can traverse the clump hierarchy to get a list of all of the 'leaf' clumps
leaf_clumps = get_lowest_clumps(master_clump)

# If you'd like to visualize these clumps, a list of clumps can be supplied to
# the "clumps" callback on a plot.  First, we create a projection plot:
prj = yt.ProjectionPlot(ds, 2, field, center='c', width=(20,'kpc'))

# Next we annotate the plot with contours on the borders of the clumps
prj.annotate_clumps(leaf_clumps)

# Lastly, we write the plot to disk.
prj.save('clumps')

# set up our namespace
from yt.mods import * 
from yt.analysis_modules.level_sets.api import *

fn = "IsolatedGalaxy/galaxy0030/galaxy0030" # parameter file to load
field = "density" # this is the field we look for contours over -- we could do
                  # this over anything.  Other common choices are 'AveragedDensity'
                  # and 'Dark_Matter_Density'.
step = 2.0 # This is the multiplicative interval between contours.

pf = load(fn) # load data

# We want to find clumps over the entire dataset, so we'll just grab the whole
# thing!  This is a convenience parameter that prepares an object that covers
# the whole domain.  Note, though, that it will load on demand and not before!
data_source = pf.disk([0.5, 0.5, 0.5], [0., 0., 1.], 
                        8./pf.units['kpc'], 1./pf.units['kpc'])

# Now we set some sane min/max values between which we want to find contours.
# This is how we tell the clump finder what to look for -- it won't look for
# contours connected below or above these threshold values.
c_min = 10**np.floor(np.log10(data_source[field]).min()  )
c_max = 10**np.floor(np.log10(data_source[field]).max()+1)

# keep only clumps with at least 20 cells
function = 'self.data[\'%s\'].size > 20' % field

# Now find get our 'base' clump -- this one just covers the whole domain.
master_clump = Clump(data_source, None, field, function=function)

# This next command accepts our base clump and we say the range between which
# we want to contour.  It recursively finds clumps within the master clump, at
# intervals defined by the step size we feed it.  The current value is
# *multiplied* by step size, rather than added to it -- so this means if you
# want to look in log10 space intervals, you would supply step = 10.0.
find_clumps(master_clump, c_min, c_max, step)

# As it goes, it appends the information about all the sub-clumps to the
# master-clump.  Among different ways we can examine it, there's a convenience
# function for outputting the full index to a file.
f = open('%s_clump_index.txt' % pf,'w')
amods.level_sets.write_clump_index(master_clump,0,f)
f.close()

# We can also output some handy information, as well.
f = open('%s_clumps.txt' % pf,'w')
amods.level_sets.write_clumps(master_clump,0,f)
f.close()

# We can traverse the clump index to get a list of all of the 'leaf' clumps
leaf_clumps = get_lowest_clumps(master_clump)

# If you'd like to visualize these clumps, a list of clumps can be supplied to
# the "clumps" callback on a plot.  First, we create a projection plot:
prj = ProjectionPlot(pf, 2, field, center='c', width=(20,'kpc'))

# Next we annotate the plot with contours on the borders of the clumps
prj.annotate_clumps(leaf_clumps)

# Lastly, we write the plot to disk.
prj.save('clumps')

# We can also save the clump object to disk to read in later so we don't have
# to spend a lot of time regenerating the clump objects.
pf.h.save_object(master_clump, 'My_clumps')

# Later, we can read in the clump object like so,
master_clump = pf.h.load_object('My_clumps')

# This script demonstrates some of the halo merger tracking infrastructure,
# for tracking halos across multiple datadumps in a time series.
# Ultimately, it outputs an HDF5 file with the important quantities for the
# top 20 most massive halos in the datadump, and it plots out their mass
# accretion histories to a series of plots.

# Currently this has only been tested with enzo outputs, but we are looking
# to generalize this with other codes imminently.

from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
from yt.analysis_modules.halo_merger_tree.api import *

# Makes a TimeSeries object from all of whatever files you have
ts = DatasetSeries.from_filenames("enzo_tiny_cosmology/DD????/DD????")

# For each datadump in our timeseries, run the friends of friends
# halo finder on it (this has only been tested with FOF currently).
# Output the information about the halos and the particles comprising each
# to disk.  These files will all be in the FOF subdirectory.
# This also works with an external FOF program run outside of yt,
# in which case skip this step and do that yourself.

# ------------------------------------------------------------
# DEPENDING ON THE SIZE OF YOUR FILES, THIS CAN BE A LONG STEP 
# but because we're writing them out to disk, you only have to do this once.
# ------------------------------------------------------------
for ds in ts:
    halo_list = FOFHaloFinder(ds)
    i = int(ds.basename[2:])
    halo_list.write_out("FOF/groups_%05i.txt" % i)
    halo_list.write_particle_lists("FOF/particles_%05i" % i)

# Create a merger tree object.  This object is a tuple, where the
# first part is a dictionary showing the correlation between file 
# output_number and redshift.  The second part is a dictionary
# correlating each halo with it's parent halos from the previous timestep
# (along with number and fraction of particles taken from that parent)

# ------------------------------------------------------------
# DEPENDING ON THE SIZE OF YOUR FILES, THIS CAN BE A LONG STEP 
# but because we're writing them out to disk, you only have to do this once.
# ------------------------------------------------------------

# by default, this saves the merger tree to disk in a CPickle:
# FOF/merger_tree.cpkl so you can use it later.
# Note that there are a bunch of filters you can place on this, so 
# that you're only building a merger tree for a subset of redshifts
# or data outputs.
mt = EnzoFOFMergerTree(external_FOF=False)

# If you want to just use your already generated merger_tree, 
# uncomment the next line.
# mt = EnzoFOFMergerTree(external_FOF=False, load_saved=True)

# For each of the top 20 most massive halos from the final timestep
# build its merger history.  You can then print this tree to screen, 
# or more usefully, save the lineage of that tree to disk for use 
# later.  save_halo_evolution() follows the largest progenitor
# which contributes the most to the subsequent children.
for i in range(20):
    mt.build_tree(i)
    mt.print_tree()
    mt.save_halo_evolution('halos.h5')

# For each of the top 20 most massive halos from the final timestep
# plot its evolution of two quantities.  The default is to look at 
# timestep vs mass, but you can look at center_of_mass phase-space 
# coordinates, fraction from progenitor, mass of halo, and more.
# These are spit out to .png files in the FOF directory.
for i in range(20):
    plot_halo_evolution('halos.h5', i)

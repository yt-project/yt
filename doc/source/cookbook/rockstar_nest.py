# You must run this job in parallel.
# There are several mpi flags which can be useful in order for it to work OK.
# It requires at least 3 processors in order to run because of the way in which
# rockstar divides up the work.  Make sure you have mpi4py installed as per
# http://yt-project.org/docs/dev/analyzing/parallel_computation.html#setting-up-parallel-yt

# Usage: mpirun -np <num_procs> --mca btl ^openib python this_script.py

import yt
from yt.analysis_modules.halo_analysis.halo_catalog import HaloCatalog
from yt.data_objects.particle_filters import add_particle_filter
from yt.analysis_modules.halo_finding.rockstar.api import RockstarHaloFinder
yt.enable_parallelism() # rockstar halofinding requires parallelism

# Create a dark matter particle filter
# This will be code dependent, but this function here is true for enzo

def DarkMatter(pfilter, data):
    filter = data[("all", "particle_type")] == 1 # DM = 1, Stars = 2
    return filter

add_particle_filter("dark_matter", function=DarkMatter, filtered_type='all', \
                    requires=["particle_type"])

# First, we make sure that this script is being run using mpirun with
# at least 3 processors as indicated in the comments above.
assert(yt.communication_system.communicators[-1].size >= 3)

# Load the dataset and apply dark matter filter
fn = "Enzo_64/DD0043/data0043"
ds = yt.load(fn)
ds.add_particle_filter('dark_matter')

# Determine highest resolution DM particle mass in sim by looking
# at the extrema of the dark_matter particle_mass field.
ad = ds.all_data()
min_dm_mass = ad.quantities.extrema(('dark_matter','particle_mass'))[0]

# Define a new particle filter to isolate all highest resolution DM particles
# and apply it to dataset
def MaxResDarkMatter(pfilter, data):
    return data["particle_mass"] <= 1.01 * min_dm_mass

add_particle_filter("max_res_dark_matter", function=MaxResDarkMatter, \
                    filtered_type='dark_matter', requires=["particle_mass"])
ds.add_particle_filter('max_res_dark_matter')

# If desired, we can see the total number of DM and High-res DM particles
#if yt.is_root():
#    print("Simulation has %d DM particles." %
#          ad['dark_matter','particle_type'].shape)
#    print("Simulation has %d Highest Res DM particles." %
#          ad['max_res_dark_matter', 'particle_type'].shape)

# Run the halo catalog on the dataset only on the highest resolution dark matter
# particles
hc = HaloCatalog(data_ds=ds, finder_method='rockstar', \
                 finder_kwargs={'dm_only':True, 'particle_type':'max_res_dark_matter'})
hc.create()

# Or alternatively, just run the RockstarHaloFinder and later import the
# output file as necessary.  You can skip this step if you've already run it
# once, but be careful since subsequent halo finds will overwrite this data.
#rhf = RockstarHaloFinder(ds, particle_type="max_res_dark_matter")
#rhf.run()
# Load the halo list from a rockstar output for this dataset
# Create a projection with the halos overplot on top
#halos = yt.load('rockstar_halos/halos_0.0.bin')
#hc = HaloCatalog(halos_ds=halos)
#hc.load()

# Regardless of your method of creating the halo catalog, use it to overplot the
# halos on a projection.
p = yt.ProjectionPlot(ds, "x", "density")
p.annotate_halos(hc, annotate_field = 'particle_identifier', width=(10,'Mpc'), factor=2)
p.save()

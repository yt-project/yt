"""
Title: Halo Mass Info
Description: This recipe finds halos and then prints out information about
             them.  Note that this recipe will take advantage of multiple CPUs
             if executed with mpirun and supplied the --parallel command line
             argument.  
Outputs: [RedshiftOutput0006_halo_info.txt]
"""
from yt.mods import *

fn = "Enzo_64/RD0006/RedshiftOutput0006" # parameter file to load
pf = load(fn) # load data

# First we run our halo finder to identify all the halos in the dataset.  This
# can take arguments, but the default are pretty sane.
halos = HaloFinder(pf)

f = open("%s_halo_info.txt" % pf, "w")

# Now, for every halo, we get the baryon data and examine it.
for halo in halos:
    # The halo has a property called 'get_sphere' that obtains a sphere
    # centered on the point of maximum density (or the center of mass, if that
    # argument is supplied) and with the radius the maximum particle radius of
    # that halo.
    sphere = halo.get_sphere()
    # We use the quantities[] method to get the total mass in baryons and in
    # particles.
    baryon_mass, particle_mass = sphere.quantities["TotalQuantity"](
            ["cell_mass", "particle_mass"])
    # Now we print out this information, along with the ID.
    f.write("Total mass in HOP group %s is %0.5e (gas = %0.5e / particles = %0.5e)\n" % \
            (halo.id, baryon_mass + particle_mass, baryon_mass, particle_mass))
f.close()

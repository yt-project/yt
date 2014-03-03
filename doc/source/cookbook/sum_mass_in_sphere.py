from yt.mods import *

# Load the dataset.
pf = load("Enzo_64/DD0029/data0029")

# Create a 1 Mpc radius sphere, centered on the max density.
sp = pf.h.sphere("max", (1.0, "mpc"))

# Use the TotalQuantity derived quantity to sum up the
# values of the cell_mass and ParticleMassMsun fields
# within the sphere.
baryon_mass, particle_mass = sp.quantities["TotalQuantity"](
        ["cell_mass", "ParticleMassMsun"])

print "Total mass in sphere is %0.5e (gas = %0.5e / particles = %0.5e)" % \
            (baryon_mass + particle_mass, baryon_mass, particle_mass)

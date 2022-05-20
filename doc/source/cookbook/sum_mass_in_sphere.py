import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0029/data0029")

# Create a 1 Mpc radius sphere, centered on the max density.
sp = ds.sphere("max", (1.0, "Mpc"))

# Use the total_quantity derived quantity to sum up the
# values of the mass and particle_mass fields
# within the sphere.
baryon_mass, particle_mass = sp.quantities.total_quantity(
    [("gas", "mass"), ("all", "particle_mass")]
)

print(
    "Total mass in sphere is %0.3e Msun (gas = %0.3e Msun, particles = %0.3e Msun)"
    % (
        (baryon_mass + particle_mass).in_units("Msun"),
        baryon_mass.in_units("Msun"),
        particle_mass.in_units("Msun"),
    )
)

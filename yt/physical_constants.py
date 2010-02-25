#
# Physical Constants and Units Conversion Factors
#

# Masses
mass_hydrogen_cgs = 1.67e-24 # g
mass_electron_cgs = 9.11e-28 # g
print "poop"
# Velocities
speed_of_light_cgs = 2.99792458e10 # cm/s, exact

# Cross Sections
cross_section_thompson_cgs = 6.65e-25 # cm^2

# Physical Constants
boltzmann_constant_cgs = 1.3806504e-16 # erg K^-1
gravitational_constant_cgs  = 6.67428e-8 # cm^3 g^-1 s^-2

rho_crit_now = 1.8788e-29 # g times h^2 (critical mass for closure, Cosmology)

# Misc. Approximations
mass_mean_atomic_cosmology = 1.22
mass_mean_atomic_galactic = 2.3

# Conversion Factors
# length
mpc_in_mpc = 1
mpc_in_kpc = 1e3
mpc_in_pc  = 1e6
mpc_in_au  = 2.063e11
mpc_in_rsun =4.43664e13
mpc_in_miles =1.917e19
mpc_in_cm  = 3.0857e24


#Short cuts
G = gravitational_constant_cgs
me = mass_electron_cgs
mp = mass_hydrogen_cgs
clight = speed_of_light_cgs
kboltz = boltzmann_constant_cgs
sigma_thompson = cross_section_thompson_cgs

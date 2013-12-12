from yt.utilities.physical_ratios import *
from yt.data_objects.yt_array import YTQuantity

mass_electron_cgs = YTQuantity(mass_electron_grams, 'g')
amu_cgs           = YTQuantity(amu_grams, 'g')
mass_hydrogen_cgs = 1.007947*amu_cgs

# Velocities
speed_of_light_cgs = YTQuantity(2.99792458e10, 'cm/s')

# Cross Sections
# 8*pi/3 (alpha*hbar*c/(2*pi))**2
cross_section_thompson_cgs = YTQuantity(6.65245854533e-25, 'cm**2')

# Charge
charge_proton_cgs = YTQuantity(4.8032056e-10, 'esu')

# Physical Constants
boltzmann_constant_cgs = YTQuantity(boltzmann_constant_erg_per_K, 'erg/K')
gravitational_constant_cgs  = YTQuantity(6.67384e-8, 'cm**3/g/s**2')
planck_constant_cgs   = YTQuantity(6.62606957e-27, 'erg*s')
stefan_boltzmann_constant_cgs = YTQuantity(5.670373e-5, 'erg/cm**2/s**1/K**4')
Tcmb = YTQuantity(2.726, 'K') # Current CMB temperature
rho_crit_now = YTQuantity(rho_crit_g_cm3_h2, 'g/cm**3/h**2') # (Cosmological critical density)

mass_sun_cgs = YTQuantity(mass_sun_grams, 'g')

#Short cuts
G = gravitational_constant_cgs
me = mass_electron_cgs
mp = mass_hydrogen_cgs
qp = charge_proton_cgs
mh = mp
clight = speed_of_light_cgs
kboltz = boltzmann_constant_cgs
hcgs = planck_constant_cgs
sigma_thompson = cross_section_thompson_cgs
Na = 1 / amu_cgs

from yt.utilities.physical_ratios import *
from yt.units.yt_array import YTQuantity
from math import pi

mass_electron_cgs = YTQuantity(mass_electron_grams, 'g')
amu_cgs           = YTQuantity(amu_grams, 'g')
mass_hydrogen_cgs = 1.007947*amu_cgs

# Velocities
speed_of_light_cgs = YTQuantity(speed_of_light_cm_per_s, 'cm/s')

# Cross Sections
# 8*pi/3 (alpha*hbar*c/(2*pi))**2
cross_section_thompson_cgs = YTQuantity(6.65245854533e-25, 'cm**2')

# Charge
charge_proton_cgs = YTQuantity(4.8032056e-10, 'esu')
elementary_charge = charge_proton_cgs

# Physical Constants
boltzmann_constant_cgs = YTQuantity(boltzmann_constant_erg_per_K, 'erg/K')
gravitational_constant_cgs  = YTQuantity(6.67384e-8, 'cm**3/g/s**2')
planck_constant_cgs   = YTQuantity(6.62606957e-27, 'erg*s')
stefan_boltzmann_constant_cgs = YTQuantity(5.670373e-5, 'erg/cm**2/s**1/K**4')
Tcmb = YTQuantity(2.726, 'K') # Current CMB temperature

# Solar System
mass_sun_cgs = YTQuantity(mass_sun_grams, 'g')
# Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
# in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
# Kluwer Academic Publishers, Dordrecht.
# REMARK: following masses include whole systems (planet + moons)
mass_jupiter_cgs = YTQuantity(mass_jupiter_grams, 'g')
mass_mercury_cgs = YTQuantity(mass_mercury_grams, 'g')
mass_venus_cgs = YTQuantity(mass_venus_grams, 'g')
mass_earth_cgs = YTQuantity(mass_earth_grams, 'g')
mass_mars_cgs = YTQuantity(mass_mars_grams, 'g')
mass_saturn_cgs = YTQuantity(mass_saturn_grams, 'g')
mass_uranus_cgs = YTQuantity(mass_uranus_grams, 'g')
mass_neptune_cgs = YTQuantity(mass_neptune_grams, 'g')

#Short cuts
G = gravitational_constant_cgs
me = mass_electron_cgs
mp = mass_hydrogen_cgs
qp = charge_proton_cgs
mh = mp
clight = speed_of_light_cgs
speed_of_light = speed_of_light_cgs
kboltz = boltzmann_constant_cgs
kb = kboltz
hcgs = planck_constant_cgs
hbar = 0.5*hcgs/pi
sigma_thompson = cross_section_thompson_cgs
Na = 1 / amu_cgs

#Planck units
m_pl = planck_mass = YTQuantity(planck_mass_grams, "g")
l_pl = planck_length = YTQuantity(planck_length_cm, "cm")
t_pl = planck_time = YTQuantity(planck_time_s, "s")
E_pl = planck_energy = YTQuantity(planck_energy_erg, "erg")
q_pl = planck_charge = YTQuantity(planck_charge_esu, "esu")
T_pl = planck_temperature = YTQuantity(planck_temperature_K, "K")

mu_0 = YTQuantity(4.0e-7*pi, "N/A**2")
eps_0 = (1.0/(clight.in_mks()**2*mu_0)).in_units("C**2/N/m**2")

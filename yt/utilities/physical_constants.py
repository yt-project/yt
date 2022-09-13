from math import pi

from yt.units.yt_array import YTQuantity
from yt.utilities.physical_ratios import (
    amu_grams,
    boltzmann_constant_erg_per_K,
    mass_earth_grams,
    mass_electron_grams,
    mass_hydrogen_grams,
    mass_jupiter_grams,
    mass_mars_grams,
    mass_mercury_grams,
    mass_neptune_grams,
    mass_saturn_grams,
    mass_sun_grams,
    mass_uranus_grams,
    mass_venus_grams,
    newton_cgs,
    planck_cgs,
    planck_charge_esu,
    planck_energy_erg,
    planck_length_cm,
    planck_mass_grams,
    planck_temperature_K,
    planck_time_s,
    speed_of_light_cm_per_s,
    standard_gravity_cm_per_s2,
)

mass_electron_cgs = YTQuantity(mass_electron_grams, "g")
mass_electron = mass_electron_cgs
me = mass_electron_cgs

amu_cgs = YTQuantity(amu_grams, "g")
amu = amu_cgs
Na = 1 / amu_cgs

mass_hydrogen_cgs = YTQuantity(mass_hydrogen_grams, "g")
mass_hydrogen = mass_hydrogen_cgs
mp = mass_hydrogen_cgs
mh = mp

# Velocities
speed_of_light_cgs = YTQuantity(speed_of_light_cm_per_s, "cm/s")
speed_of_light = speed_of_light_cgs
clight = speed_of_light_cgs
c = speed_of_light_cgs

# Cross Sections
# 8*pi/3 (alpha*hbar*c/(2*pi))**2
cross_section_thompson_cgs = YTQuantity(6.65245854533e-25, "cm**2")
cross_section_thompson = cross_section_thompson_cgs
thompson_cross_section = cross_section_thompson_cgs
sigma_thompson = cross_section_thompson_cgs

# Charge
charge_proton_cgs = YTQuantity(4.8032056e-10, "esu")
charge_proton = charge_proton_cgs
proton_charge = charge_proton_cgs
elementary_charge = charge_proton_cgs
qp = charge_proton_cgs

# Physical Constants
boltzmann_constant_cgs = YTQuantity(boltzmann_constant_erg_per_K, "erg/K")
boltzmann_constant = boltzmann_constant_cgs
kboltz = boltzmann_constant_cgs
kb = kboltz

gravitational_constant_cgs = YTQuantity(newton_cgs, "cm**3/g/s**2")
gravitational_constant = gravitational_constant_cgs
G = gravitational_constant_cgs

planck_constant_cgs = YTQuantity(planck_cgs, "erg*s")
planck_constant = planck_constant_cgs
hcgs = planck_constant_cgs
hbar = 0.5 * hcgs / pi

stefan_boltzmann_constant_cgs = YTQuantity(5.670373e-5, "erg/cm**2/s**1/K**4")
stefan_boltzmann_constant = stefan_boltzmann_constant_cgs

Tcmb = YTQuantity(2.726, "K")  # Current CMB temperature
CMB_temperature = Tcmb

# Solar System
mass_sun_cgs = YTQuantity(mass_sun_grams, "g")
mass_sun = mass_sun_cgs
solar_mass = mass_sun_cgs
msun = mass_sun_cgs
# Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards",
# in Highlights of Astronomy (I. Appenzeller, ed.), Table 1,
# Kluwer Academic Publishers, Dordrecht.
# REMARK: following masses include whole systems (planet + moons)
mass_jupiter_cgs = YTQuantity(mass_jupiter_grams, "g")
mass_jupiter = mass_jupiter_cgs
jupiter_mas = mass_jupiter_cgs

mass_mercury_cgs = YTQuantity(mass_mercury_grams, "g")
mass_mercury = mass_mercury_cgs
mercury_mass = mass_mercury_cgs

mass_venus_cgs = YTQuantity(mass_venus_grams, "g")
mass_venus = mass_venus_cgs
venus_mass = mass_venus_cgs

mass_earth_cgs = YTQuantity(mass_earth_grams, "g")
mass_earth = mass_earth_cgs
earth_mass = mass_earth_cgs
mearth = mass_earth_cgs

mass_mars_cgs = YTQuantity(mass_mars_grams, "g")
mass_mars = mass_mars_cgs
mars_mass = mass_mars_cgs

mass_saturn_cgs = YTQuantity(mass_saturn_grams, "g")
mass_saturn = mass_saturn_cgs
saturn_mass = mass_saturn_cgs

mass_uranus_cgs = YTQuantity(mass_uranus_grams, "g")
mass_uranus = mass_uranus_cgs
uranus_mass = mass_uranus_cgs

mass_neptune_cgs = YTQuantity(mass_neptune_grams, "g")
mass_neptune = mass_neptune_cgs
neptune_mass = mass_neptune_cgs

# Planck units
m_pl = planck_mass = YTQuantity(planck_mass_grams, "g")
l_pl = planck_length = YTQuantity(planck_length_cm, "cm")
t_pl = planck_time = YTQuantity(planck_time_s, "s")
E_pl = planck_energy = YTQuantity(planck_energy_erg, "erg")
q_pl = planck_charge = YTQuantity(planck_charge_esu, "esu")
T_pl = planck_temperature = YTQuantity(planck_temperature_K, "K")

# MKS E&M units
mu_0 = YTQuantity(4.0e-7 * pi, "N/A**2")
eps_0 = (1.0 / (clight**2 * mu_0)).in_units("C**2/N/m**2")

# Misc
standard_gravity_cgs = YTQuantity(standard_gravity_cm_per_s2, "cm/s**2")
standard_gravity = standard_gravity_cgs

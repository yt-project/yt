"""
The default unit symbol lookup table.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.units import dimensions
from yt.utilities.physical_ratios import \
    cm_per_pc, cm_per_ly, cm_per_au, cm_per_rsun, \
    mass_sun_grams, sec_per_year, sec_per_day, sec_per_hr, \
    sec_per_min, temp_sun_kelvin, luminosity_sun_ergs_per_sec, \
    metallicity_sun, erg_per_eV, amu_grams, mass_electron_grams, \
    cm_per_ang, jansky_cgs, mass_jupiter_grams, mass_earth_grams, \
    boltzmann_constant_erg_per_K, kelvin_per_rankine, \
    speed_of_light_cm_per_s, planck_length_cm, planck_charge_esu, \
    planck_energy_erg, planck_mass_grams, planck_temperature_K, planck_time_s
import numpy as np

# Lookup a unit symbol with the symbol string, and provide a tuple with the
# conversion factor to cgs and dimensionality.

default_unit_symbol_lut = {
    # base
    "g":  (1.0, dimensions.mass),
    #"cm": (1.0, length, r"\rm{cm}"),  # duplicate with meter below...
    "s":  (1.0, dimensions.time),
    "K":  (1.0, dimensions.temperature),
    "radian": (1.0, dimensions.angle),

    # other cgs
    "dyne": (1.0, dimensions.force),
    "erg":  (1.0, dimensions.energy),
    "esu":  (1.0, dimensions.charge_cgs),
    "gauss": (1.0, dimensions.magnetic_field_cgs),
    "degC": (1.0, dimensions.temperature, -273.15),
    "statA": (1.0, dimensions.current_cgs),

    # some SI
    "m": (1.0e2, dimensions.length),
    "J": (1.0e7, dimensions.energy),
    "W": (1.0e7, dimensions.power),
    "Hz": (1.0, dimensions.rate),
    "N": (1.0e5, dimensions.force),
    "C": (0.1*speed_of_light_cm_per_s, dimensions.charge_mks),
    "A": (0.1*speed_of_light_cm_per_s, dimensions.current_mks),
    "T": (1.0e4, dimensions.magnetic_field_mks),
    "Pa": (10.0, dimensions.pressure),

    # Imperial units
    "ft": (30.48, dimensions.length),
    "mile": (160934, dimensions.length),
    "degF": (kelvin_per_rankine, dimensions.temperature, -459.67),
    "R": (kelvin_per_rankine, dimensions.temperature),

    # dimensionless stuff
    "h": (1.0, dimensions.dimensionless), # needs to be added for rho_crit_now
    "dimensionless": (1.0, dimensions.dimensionless),

    # times
    "min": (sec_per_min, dimensions.time),
    "hr":  (sec_per_hr, dimensions.time),
    "day": (sec_per_day, dimensions.time),
    "yr":  (sec_per_year, dimensions.time),

    # Velocities
    "c": (speed_of_light_cm_per_s, dimensions.velocity),

    # Solar units
    "Msun": (mass_sun_grams, dimensions.mass),
    "msun": (mass_sun_grams, dimensions.mass),
    "Rsun": (cm_per_rsun, dimensions.length),
    "rsun": (cm_per_rsun, dimensions.length),
    "Lsun": (luminosity_sun_ergs_per_sec, dimensions.power),
    "Tsun": (temp_sun_kelvin, dimensions.temperature),
    "Zsun": (metallicity_sun, dimensions.dimensionless),
    "Mjup": (mass_jupiter_grams, dimensions.mass),
    "Mearth": (mass_earth_grams, dimensions.mass),

    # astro distances
    "AU": (cm_per_au, dimensions.length),
    "au": (cm_per_au, dimensions.length),
    "ly": (cm_per_ly, dimensions.length),
    "pc": (cm_per_pc, dimensions.length),

    # angles
    "degree": (np.pi/180., dimensions.angle), # degrees
    "arcmin": (np.pi/10800., dimensions.angle), # arcminutes
    "arcsec": (np.pi/648000., dimensions.angle), # arcseconds
    "mas": (np.pi/648000000., dimensions.angle), # millarcseconds
    "hourangle": (np.pi/12., dimensions.angle), # hour angle
    "steradian": (1.0, dimensions.solid_angle),

    # misc
    "eV": (erg_per_eV, dimensions.energy),
    "amu": (amu_grams, dimensions.mass),
    "angstrom": (cm_per_ang, dimensions.length),
    "Jy": (jansky_cgs, dimensions.specific_flux),
    "counts": (1.0, dimensions.dimensionless),
    "photons": (1.0, dimensions.dimensionless),
    "me": (mass_electron_grams, dimensions.mass),

    # for AstroPy compatibility
    "solMass": (mass_sun_grams, dimensions.mass),
    "solRad": (cm_per_rsun, dimensions.length),
    "solLum": (luminosity_sun_ergs_per_sec, dimensions.power),
    "dyn": (1.0, dimensions.force),
    "sr": (1.0, dimensions.solid_angle),
    "rad": (1.0, dimensions.solid_angle),
    "deg": (np.pi/180., dimensions.angle),
    "Fr":  (1.0, dimensions.charge_cgs),
    "G": (1.0, dimensions.magnetic_field_cgs),
    "d": (1.0, dimensions.time),
    "Angstrom": (cm_per_ang, dimensions.length),

    # Planck units
    "m_pl": (planck_mass_grams, dimensions.mass),
    "l_pl": (planck_length_cm, dimensions.length),
    "t_pl": (planck_time_s, dimensions.time),
    "T_pl": (planck_temperature_K, dimensions.temperature),
    "q_pl": (planck_charge_esu, dimensions.charge_cgs),
    "E_pl": (planck_energy_erg, dimensions.energy),

}

# Add LaTeX representations for units with trivial representations.
latex_symbol_lut = {
    "unitary" : r"",
    "dimensionless" : r"",
    "code_length" : r"\rm{code}\ \rm{length}",
    "code_time" : r"\rm{code}\ \rm{time}",
    "code_mass" : r"\rm{code}\ \rm{mass}",
    "code_temperature" : r"\rm{code}\ \rm{temperature}",
    "code_metallicity" : r"\rm{code}\ \rm{metallicity}",
    "code_velocity" : r"\rm{code}\ \rm{velocity}",
    "code_magnetic" : r"\rm{code}\ \rm{magnetic}",
    "Msun" : r"\rm{M}_\odot",
    "msun" : r"\rm{M}_\odot",
    "Rsun" : r"\rm{R}_\odot",
    "rsun" : r"\rm{R}_\odot",
    "Lsun" : r"\rm{L}_\odot",
    "Tsun" : r"\rm{T}_\odot",
    "Zsun" : r"\rm{Z}_\odot",
}
for key in default_unit_symbol_lut:
    if key not in latex_symbol_lut:
        latex_symbol_lut[key] = r"\rm{" + key + r"}"

# This dictionary formatting from magnitude package, credit to Juan Reyero.
unit_prefixes = {
    'Y': 1e24,   # yotta
    'Z': 1e21,   # zetta
    'E': 1e18,   # exa
    'P': 1e15,   # peta
    'T': 1e12,   # tera
    'G': 1e9,    # giga
    'M': 1e6,    # mega
    'k': 1e3,    # kilo
    'd': 1e1,    # deci
    'c': 1e-2,   # centi
    'm': 1e-3,   # mili
    'u': 1e-6,   # micro
    'n': 1e-9,   # nano
    'p': 1e-12,  # pico
    'f': 1e-15,  # femto
    'a': 1e-18,  # atto
    'z': 1e-21,  # zepto
    'y': 1e-24,  # yocto
}

latex_prefixes = {
    "u" : "\\mu",
    }

prefixable_units = (
    "m",
    "pc",
    "mcm",
    "pccm",
    "g",
    "eV",
    "s",
    "yr",
    "K",
    "dyne",
    "erg",
    "esu",
    "J",
    "Hz",
    "W",
    "gauss",
    "G",
    "Jy",
    "N",
    "T",
    "A",
    "C",
    "statA",
    "Pa",
)

cgs_base_units = {
    dimensions.mass:'g',
    dimensions.length:'cm',
    dimensions.time:'s',
    dimensions.temperature:'K',
    dimensions.angle:'radian',
}

mks_base_units = {
    dimensions.mass:'kg',
    dimensions.length:'m',
    dimensions.time:'s',
    dimensions.temperature:'K',
    dimensions.angle:'radian',
    dimensions.current_mks:'A',
}

cgs_conversions = {
    "C":"esu",
    "T":"gauss",
    "A":"statA",
}

for conv in cgs_conversions.keys():
    if conv in prefixable_units:
        for p in unit_prefixes:
            cgs_conversions[p+conv] = p+cgs_conversions[conv]

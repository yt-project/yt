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
    hubble_constant_hertz
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

    # "code" units, default to CGS conversion.
    # These default values are overridden in the code frontends
    "code_length" : (1.0, dimensions.length),
    "unitary"   : (1.0, dimensions.length),
    "code_mass" : (1.0, dimensions.mass),
    "code_time" : (1.0, dimensions.time),
    "code_velocity" : (1.0, dimensions.velocity),
    "code_magnetic" : (1.0, dimensions.magnetic_field),
    "code_temperature" : (1.0, dimensions.temperature),
    "code_metallicity" : (1.0, dimensions.metallicity),

    # other cgs
    "dyne": (1.0, dimensions.force),
    "erg":  (1.0, dimensions.energy),
    "esu":  (1.0, dimensions.charge),
    "gauss": (1.0, dimensions.magnetic_field),

    # some SI
    "m": (1.0e2, dimensions.length),
    "J": (1.0e7, dimensions.energy),
    "W": (1.0e7, dimensions.power),
    "Hz": (1.0, dimensions.rate),

    # Imperial units
    "ft": (30.48, dimensions.length),
    "mile": (160934, dimensions.length),

    # dimensionless stuff
    "h": (1.0, dimensions.dimensionless), # needs to be added for rho_crit_now

    # times
    "min": (sec_per_min, dimensions.time),
    "hr":  (sec_per_hr, dimensions.time),
    "day": (sec_per_day, dimensions.time),
    "yr":  (sec_per_year, dimensions.time),

    # Solar units
    "Msun": ( mass_sun_grams, dimensions.mass),
    "msun": ( mass_sun_grams, dimensions.mass),
    "Rsun": ( cm_per_rsun, dimensions.length),
    "rsun": ( cm_per_rsun, dimensions.length),
    "Lsun": ( luminosity_sun_ergs_per_sec, dimensions.power),
    "Tsun": ( temp_sun_kelvin, dimensions.temperature),
    "Zsun": ( metallicity_sun, dimensions.metallicity),

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
    
    # misc
    "eV": (erg_per_eV, dimensions.energy),
    "amu": (amu_grams, dimensions.mass),
    "me": (mass_electron_grams, dimensions.mass),

}

# Add LaTeX representations for units with trivial representations.
latex_symbol_lut = {
    "unitary" : "",
    "code_length" : "\\rm{code}\/\\rm{length}",
    "code_time" : "\\rm{code}\/\\rm{time}",
    "code_mass" : "\\rm{code}\/\\rm{mass}",
    "code_temperature" : "\\rm{code}\/\\rm{temperature}",
    "code_metallicity" : "\\rm{code}\/\\rm{metallicity}",
    "code_velocity" : "\\rm{code}\/\\rm{velocity}",
    "Msun" : "\\rm{M}_\\odot",
    "msun" : "\\rm{M}_\\odot",
    "Rsun" : "\\rm{R}_\\odot",
    "rsun" : "\\rm{R}_\\odot",
    "Lsun" : "\\rm{L}_\\odot",
    "Tsun" : "\\rm{T}_\\odot",
    "Zsun" : "\\rm{Z}_\\odot",
}
for key in default_unit_symbol_lut:
    if key not in latex_symbol_lut:
        latex_symbol_lut[key] = "\\rm{" + key + "}"

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
)

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
    cm_per_pc, cm_per_ly, cm_per_au, cm_per_rsun, cm_per_m, \
    mass_sun_grams, sec_per_year, sec_per_day, sec_per_hr, \
    sec_per_min, temp_sun_kelvin, luminosity_sun_ergs_per_sec, \
    metallicity_sun, erg_per_eV, amu_grams, mass_electron_grams, \
    cm_per_ang, jansky_cgs, mass_jupiter_grams, mass_earth_grams, \
    kelvin_per_rankine, speed_of_light_cm_per_s, planck_length_cm, \
    planck_charge_esu, planck_energy_erg, planck_mass_grams, \
    planck_temperature_K, planck_time_s, mass_hydrogen_grams, \
    grams_per_pound, standard_gravity_cm_per_s2, pascal_per_atm, \
    newton_cgs, cm_per_rearth, cm_per_rjup
import numpy as np

# Lookup a unit symbol with the symbol string, and provide a tuple with the
# conversion factor to cgs and dimensionality.

default_unit_symbol_lut = {
    # base
    "g":  (1.0, dimensions.mass, 0.0, r"\rm{g}"),
    "s":  (1.0, dimensions.time, 0.0, r"\rm{s}"),
    "K":  (1.0, dimensions.temperature, 0.0, r"\rm{K}"),
    "radian": (1.0, dimensions.angle, 0.0, r"\rm{radian}"),

    # other cgs
    "dyne": (1.0, dimensions.force, 0.0, r"\rm{dyn}"),
    "erg":  (1.0, dimensions.energy, 0.0, r"\rm{erg}"),
    "esu":  (1.0, dimensions.charge_cgs, 0.0, r"\rm{esu}"),
    "gauss": (1.0, dimensions.magnetic_field_cgs, 0.0, r"\rm{G}"),
    "degC": (1.0, dimensions.temperature, -273.15, r"^\circ\rm{C}"),
    "statA": (1.0, dimensions.current_cgs, 0.0, r"\rm{statA}"),
    "statV": (1.0, dimensions.electric_potential_cgs, 0.0, r"\rm{statV}"),
    "statohm": (1.0, dimensions.resistance_cgs, 0.0, r"\rm{statohm}"),

    # some SI
    "m": (1.0e2, dimensions.length, 0.0, r"\rm{m}"),
    "J": (1.0e7, dimensions.energy, 0.0, r"\rm{J}"),
    "W": (1.0e7, dimensions.power, 0.0, r"\rm{W}"),
    "Hz": (1.0, dimensions.rate, 0.0, r"\rm{Hz}"),
    "N": (1.0e5, dimensions.force, 0.0, r"\rm{N}"),
    "C": (1.0, dimensions.charge_mks, 0.0, r"\rm{C}"),
    "A": (1.0, dimensions.current_mks, 0.0, r"\rm{A}"),
    "T": (1000.0, dimensions.magnetic_field_mks, 0.0, r"\rm{T}"),
    "Pa": (10.0, dimensions.pressure, 0.0, r"\rm{Pa}"),
    "V": (1.0e7, dimensions.electric_potential_mks, 0.0, r"\rm{V}"),
    "ohm": (1.0e7, dimensions.resistance_mks, 0.0, r"\Omega"),

    # Imperial and other non-metric units
    "ft": (30.48, dimensions.length, 0.0, r"\rm{ft}"),
    "mile": (160934, dimensions.length, 0.0, r"\rm{mile}"),
    "degF": (kelvin_per_rankine, dimensions.temperature, -459.67,
             "^\circ\rm{F}"),
    "R": (kelvin_per_rankine, dimensions.temperature, 0.0, r"^\circ\rm{R}"),
    "lbf": (grams_per_pound*standard_gravity_cm_per_s2, dimensions.force, 0.0, r"\rm{lbf}"),
    "lbm": (grams_per_pound, dimensions.mass, 0.0, r"\rm{lbm}"),
    "atm": (pascal_per_atm*10., dimensions.pressure, 0.0, r"\rm{atm}"),

    # dimensionless stuff
    "h": (1.0, dimensions.dimensionless, 0.0, r"h"),  # needs to be added for rho_crit_now
    "dimensionless": (1.0, dimensions.dimensionless, 0.0, r""),

    # times
    "min": (sec_per_min, dimensions.time, 0.0, r"\rm{min}"),
    "hr":  (sec_per_hr, dimensions.time, 0.0, r"\rm{hr}"),
    "day": (sec_per_day, dimensions.time, 0.0, r"\rm{d}"),
    "yr":  (sec_per_year, dimensions.time, 0.0, r"\rm{yr}"),

    # Velocities
    "c": (speed_of_light_cm_per_s, dimensions.velocity, 0.0, r"\rm{c}"),

    # Solar units
    "Msun": (mass_sun_grams, dimensions.mass, 0.0, r"M_\odot"),
    "msun": (mass_sun_grams, dimensions.mass, 0.0, r"M_\odot"),
    "Rsun": (cm_per_rsun, dimensions.length, 0.0, r"R_\odot"),
    "rsun": (cm_per_rsun, dimensions.length, 0.0, r"R_\odot"),
    "R_sun": (cm_per_rsun, dimensions.length, 0.0, r"R_\odot"),
    "r_sun": (cm_per_rsun, dimensions.length, 0.0, r"R_\odot"),
    "Lsun": (luminosity_sun_ergs_per_sec, dimensions.power, 0.0, r"L_\odot"),
    "Tsun": (temp_sun_kelvin, dimensions.temperature, 0.0, r"T_\odot"),
    "Zsun": (metallicity_sun, dimensions.dimensionless, 0.0, r"Z_\odot"),
    "Mjup": (mass_jupiter_grams, dimensions.mass, 0.0, r"M_{\rm{Jup}}"),
    "Mearth": (mass_earth_grams, dimensions.mass, 0.0, r"M_\oplus"),

    # astro distances
    "AU": (cm_per_au, dimensions.length, 0.0, r"\rm{AU}"),
    "au": (cm_per_au, dimensions.length, 0.0, r"\rm{AU}"),
    "ly": (cm_per_ly, dimensions.length, 0.0, r"\rm{ly}"),
    "pc": (cm_per_pc, dimensions.length, 0.0, r"\rm{pc}"),

    # angles
    "degree": (np.pi/180., dimensions.angle, 0.0, r"\rm{deg}"),  # degrees
    "arcmin": (np.pi/10800., dimensions.angle, 0.0,
               r"\rm{arcmin}"),  # arcminutes
    "arcsec": (np.pi/648000., dimensions.angle, 0.0,
               r"\rm{arcsec}"),  # arcseconds
    "mas": (np.pi/648000000., dimensions.angle, 0.0,
            r"\rm{mas}"),  # milliarcseconds
    "hourangle": (np.pi/12., dimensions.angle, 0.0, r"\rm{HA}"),  # hour angle
    "steradian": (1.0, dimensions.solid_angle, 0.0, r"\rm{sr}"),
    "lat": (-np.pi/180.0, dimensions.angle, 90.0, r"\rm{Latitude}"),
    "lon": (np.pi/180.0, dimensions.angle, -180.0, r"\rm{Longitude}"),

    # misc
    "eV": (erg_per_eV, dimensions.energy, 0.0, r"\rm{eV}"),
    "amu": (amu_grams, dimensions.mass, 0.0, r"\rm{amu}"),
    "angstrom": (cm_per_ang, dimensions.length, 0.0, r"\AA"),
    "Jy": (jansky_cgs, dimensions.specific_flux, 0.0, r"\rm{Jy}"),
    "counts": (1.0, dimensions.dimensionless, 0.0, r"\rm{counts}"),
    "photons": (1.0, dimensions.dimensionless, 0.0, r"\rm{photons}"),
    "me": (mass_electron_grams, dimensions.mass, 0.0, r"m_e"),
    "mp": (mass_hydrogen_grams, dimensions.mass, 0.0, r"m_p"),
    "mol": (1.0/amu_grams, dimensions.dimensionless, 0.0, r"\rm{mol}"),
    'Sv': (cm_per_m**2, dimensions.specific_energy, 0.0,
           r"\rm{Sv}"),

    # for AstroPy compatibility
    "solMass": (mass_sun_grams, dimensions.mass, 0.0, r"M_\odot"),
    "solRad": (cm_per_rsun, dimensions.length, 0.0, r"R_\odot"),
    "solLum": (luminosity_sun_ergs_per_sec, dimensions.power, 0.0, r"L_\odot"),
    "dyn": (1.0, dimensions.force, 0.0, r"\rm{dyn}"),
    "sr": (1.0, dimensions.solid_angle, 0.0, r"\rm{sr}"),
    "rad": (1.0, dimensions.solid_angle, 0.0, r"\rm{rad}"),
    "deg": (np.pi/180., dimensions.angle, 0.0, r"\rm{deg}"),
    "Fr":  (1.0, dimensions.charge_cgs, 0.0, r"\rm{Fr}"),
    "G": (1.0, dimensions.magnetic_field_cgs, 0.0, r"\rm{G}"),
    "d": (1.0, dimensions.time, 0.0, r"\rm{d}"),
    "Angstrom": (cm_per_ang, dimensions.length, 0.0, r"\AA"),
    "statC": (1.0, dimensions.charge_cgs, 0.0, r"\rm{statC}"),

    # Planck units
    "m_pl": (planck_mass_grams, dimensions.mass, 0.0, r"m_{\rm{P}}"),
    "l_pl": (planck_length_cm, dimensions.length, 0.0, r"\ell_\rm{P}"),
    "t_pl": (planck_time_s, dimensions.time, 0.0, r"t_{\rm{P}}"),
    "T_pl": (planck_temperature_K, dimensions.temperature, 0.0, r"T_{\rm{P}}"),
    "q_pl": (planck_charge_esu, dimensions.charge_cgs, 0.0, r"q_{\rm{P}}"),
    "E_pl": (planck_energy_erg, dimensions.energy, 0.0, r"E_{\rm{P}}"),

    # Geometrized units
    "m_geom": (mass_sun_grams, dimensions.mass, 0.0, r"M_\odot"),
    "l_geom": (newton_cgs*mass_sun_grams/speed_of_light_cm_per_s**2, dimensions.length, 0.0, r"M_\odot"),
    "t_geom": (newton_cgs*mass_sun_grams/speed_of_light_cm_per_s**3, dimensions.time, 0.0, r"M_\odot"),

    # Some Solar System units
    "R_earth": (cm_per_rearth, dimensions.length, 0.0, r"R_\oplus"),
    "r_earth": (cm_per_rearth, dimensions.length, 0.0, r"R_\oplus"),
    "R_jup": (cm_per_rjup, dimensions.length, 0.0, r"R_\mathrm{Jup}"),
    "r_jup": (cm_per_rjup, dimensions.length, 0.0, r"R_\mathrm{Jup}"),
}

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
    "u": r"\mu",
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
    "V",
    "statV",
    "ohm",
    "statohm",
    "Sv",
)

default_base_units = {
    dimensions.mass: 'g',
    dimensions.length: 'cm',
    dimensions.time: 's',
    dimensions.temperature: 'K',
    dimensions.angle: 'radian',
    dimensions.current_mks: 'A',
}

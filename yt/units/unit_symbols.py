"""
A place to statically create unit quantities.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.units.yt_array import YTQuantity as quan

#
# meter
#

fm = femtometer = quan(1.0, "fm")
pm = picometer  = quan(1.0, "pm")
nm = nanometer  = quan(1.0, "nm")
um = micrometer = quan(1.0, "um")
mm = millimeter = quan(1.0, "mm")
cm = centimeter = quan(1.0, "cm")
m  = meter      = quan(1.0, "m")
km = kilometer  = quan(1.0, "km")
Mm = Megameter  = megameter = quan(1.0, "Mm")

#
# parsec
#

pc  = parsec = quan(1.0, "pc")
kpc = kiloparsec = quan(1.0, "kpc")
Mpc = mpc = megaparsec = quan(1.0, "Mpc")
Gpc = gpc = Gigaparsec = quan(1.0, "Gpc")

#
# gram
#

mg = milligram  = quan(1.0, "mg")
g  = gram       = quan(1.0, "g")
kg = kilogram   = quan(1.0, "kg")

#
# second
#

fs   = femtoseconds = quan(1.0, "fs")
ps   = picosecond   = quan(1.0, "ps")
ns   = nanosecond   = quan(1.0, "ns")
ms   = millisecond  = quan(1.0, "ms")
s    = second       = quan(1.0, "s")

#
# minute
#

min = minute = quan(1.0, "min")

#
# hr
#

hr = hour = quan(1.0, "hr")

#
# day
#

day = quan(1.0, "day")

#
# year
#

yr   = year                = quan(1.0, "yr")
kyr  = kiloyear            = quan(1.0, "kyr")
Myr  = Megayear = megayear = quan(1.0, "Myr")
Gyr  = Gigayear = gigayear = quan(1.0, "Gyr")

#
# Kelvin
#

degree_kelvin = Kelvin = K = quan(1.0, "K")

#
# Misc CGS
#

dyne = dyn = quan(1.0, "dyne")
erg = ergs = quan(1.0, "erg")

#
# Misc SI
#

N  = Newton = newton = quan(1.0, "N")
J  = Joule = joule = quan(1.0, "J")
W  = Watt  = watt = quan(1.0, "W")
Hz = Hertz = hertz = quan(1.0, "Hz") 

#
# Imperial units
#

ft = foot = quan(1.0, "ft")
mile = quan(1.0, "mile")

#
# Solar units
#

Msun = solar_mass = quan(1.0, "Msun")
msun = quan(1.0, "msun")
Rsun = R_sun = solar_radius = quan(1.0, "Rsun")
rsun = r_sun = quan(1.0, "rsun")
Lsun = lsun = l_sun = solar_luminosity = quan(1.0, "Lsun")
Tsun = T_sun = solar_temperature = quan(1.0, "Tsun")
Zsun = Z_sun = solar_metallicity = quan(1.0, "Zsun")

#
# Misc Astronomical units
#

AU = astronomical_unit = quan(1.0, "AU")
au = quan(1.0, "au")
ly = light_year = quan(1.0, "ly")
Rearth = R_earth = earth_radius = quan(1.0, 'R_earth')
rearth = r_earth = quan(1.0, 'r_earth')
Rjup = R_jup = jupiter_radius = quan(1.0, 'R_jup')
rjup = r_jup = quan(1.0, 'r_jup')

#
# Physical units
#

eV  = electron_volt = quan(1.0, "eV")
keV = kilo_electron_volt = quan(1.0, "keV")
MeV = mega_electron_volt = quan(1.0, "MeV")
GeV = giga_electron_volt = quan(1.0, "GeV")
amu = atomic_mass_unit = quan(1.0, "amu")
angstrom = quan(1.0, "angstrom")
me  = electron_mass = quan(1.0, "me")

#
# Angle units
#

deg    = degree = quan(1.0, "degree")
rad    = radian = quan(1.0, "radian")
arcsec = arcsecond = quan(1.0, "arcsec")
arcmin = arcminute = quan(1.0, "arcmin")
mas    = milliarcsecond = quan(1.0, "mas")
sr     = steradian = quan(1.0, "steradian")

#
# CGS electromagnetic units
#

electrostatic_unit = esu = quan(1.0, "esu")
gauss = G = quan(1.0, "gauss")
statampere = statA = quan(1.0, "statA")
statvolt = statV = quan(1.0, "statV")
statohm = quan(1.0, "statohm")

#
# SI electromagnetic units
#

C = coulomb = Coulomb = quan(1.0, "C")
T = tesla = Tesla = quan(1.0, "T")
A = ampere = Ampere = quan(1.0, "A")
V = volt = Volt = quan(1.0, "V")
ohm = Ohm = quan(1.0, "ohm")

#
# Geographic units
#

latitude = lat = quan(1.0, "lat")
longitude = lon = quan(1.0, "lon")

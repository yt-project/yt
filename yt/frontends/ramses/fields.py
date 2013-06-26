"""
RAMSES-specific fields

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    NullFunc, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
from yt.data_objects.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs, \
    mass_sun_cgs, \
    mh
import numpy as np

RAMSESFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo, "RFI")
add_field = RAMSESFieldInfo.add_field

KnownRAMSESFields = FieldInfoContainer()
add_ramses_field = KnownRAMSESFields.add_field

known_ramses_fields = [
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
]

for f in known_ramses_fields:
    if f not in KnownRAMSESFields:
        add_ramses_field(f, function=NullFunc, take_log=True,
                  validators = [ValidateDataField(f)])

def dx(field, data):
    return data.fwidth[:,0]
add_field("dx", function=dx)

def dy(field, data):
    return data.fwidth[:,1]
add_field("dy", function=dy)

def dz(field, data):
    return data.fwidth[:,2]
add_field("dz", function=dz)

def _convertDensity(data):
    return data.convert("Density")
KnownRAMSESFields["Density"].units = "g / cm**3"
KnownRAMSESFields["Density"]._convert_function=_convertDensity

def _convertPressure(data):
    return data.convert("Pressure")
KnownRAMSESFields["Pressure"]._units="dyne/cm**2"
KnownRAMSESFields["Pressure"]._convert_function=_convertPressure

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = KnownRAMSESFields["%s-velocity" % ax]
    f.units = "cm / s"
    f._convert_function = _convertVelocity
    f.take_log = False

known_ramses_particle_fields = [
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_velocity_x",
    "particle_velocity_y",
    "particle_velocity_z",
    "particle_mass",
    "particle_identifier",
    "particle_refinement_level",
    "particle_age",
    "particle_metallicity",
]

for f in known_ramses_particle_fields:
    add_ramses_field(("all", f), function=NullFunc, take_log=True,
              particle_type = True)

for ax in 'xyz':
    KnownRAMSESFields["all", "particle_velocity_%s" % ax]._convert_function = \
        _convertVelocity

def _convertParticleMass(data):
    return data.convert("mass")

KnownRAMSESFields["all", "particle_mass"]._convert_function = \
        _convertParticleMass
KnownRAMSESFields["all", "particle_mass"]._units = "g"

def _Temperature(field, data):
    rv = data["Pressure"]/data["Density"]
    rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
    return rv
add_field("Temperature", function=_Temperature, units="K")

# We'll add a bunch of species fields here.  In the not too distant future,
# we'll be moving all of these to a unified field location, so they can be
# shared between various frontends.

# NOTE: No Electron here because I don't know how RAMSES handles them, and if
# they are handled differently than Enzo does (where they are scaled to mh)

_speciesList = ["HI", "HII",
                "HeI", "HeII", "HeIII",
                "H2I", "H2II", "HM",
                "DI", "DII", "HDI"]
_speciesMass = {"HI": 1.0, "HII": 1.0,
                "HeI": 4.0, "HeII": 4.0, "HeIII": 4.0,
                "H2I": 2.0, "H2II": 2.0, "HM": 1.0,
                "DI": 2.0, "DII": 2.0, "HDI": 3.0}

def _SpeciesComovingDensity(field, data):
    sp = field.name.split("_")[0] + "_Density"
    ef = (1.0 + data.pf.current_redshift)**3.0
    return data[sp] / ef

def _SpeciesFraction(field, data):
    sp = field.name.split("_")[0] + "_Density"
    return data[sp] / data["Density"]

def _SpeciesMass(field, data):
    sp = field.name.split("_")[0] + "_Density"
    return data[sp] * data["CellVolume"]

def _SpeciesNumberDensity(field, data):
    species = field.name.split("_")[0]
    sp = field.name.split("_")[0] + "_Density"
    return data[sp] / _speciesMass[species]

def _convertCellMassMsun(data):
    return 1.0/mass_sun_cgs # g^-1
def _ConvertNumberDensity(data):
    return 1.0/mh

for species in _speciesList:
    add_ramses_field("%s_Density" % species,
             function = NullFunc,
             display_name = "%s\/Density" % species,
             convert_function = _convertDensity,
             units = "g/cm**3")
    add_field("%s_Fraction" % species,
             function=_SpeciesFraction,
             validators=ValidateDataField("%s_Density" % species),
             display_name="%s\/Fraction" % species)
    add_field("Comoving_%s_Density" % species,
             function=_SpeciesComovingDensity,
             validators=ValidateDataField("%s_Density" % species),
             display_name="Comoving\/%s\/Density" % species)
    add_field("%s_Mass" % species, units="g", 
              function=_SpeciesMass, 
              validators=ValidateDataField("%s_Density" % species),
              display_name="%s\/Mass" % species)
    add_field("%s_MassMsun" % species, units="Msun",
              function=_SpeciesMass, 
              convert_function=_convertCellMassMsun,
              validators=ValidateDataField("%s_Density" % species),
              display_name="%s\/Mass" % species)
    if _speciesMass.has_key(species):
        add_field("%s_NumberDensity" % species,
                  function=_SpeciesNumberDensity,
                  convert_function=_ConvertNumberDensity,
                  validators=ValidateDataField("%s_Density" % species))

# PARTICLE FIELDS
particle_vector_functions("all", ["particle_position_%s" % ax for ax in 'xyz'],
                                 ["particle_velocity_%s" % ax for ax in 'xyz'],
                          RAMSESFieldInfo)
particle_deposition_functions("all", "Coordinates", "particle_mass",
                               RAMSESFieldInfo)

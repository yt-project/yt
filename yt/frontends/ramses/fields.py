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
KnownRAMSESFields["Density"]._units = r"\rm{g}/\rm{cm}^3"
KnownRAMSESFields["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownRAMSESFields["Density"]._convert_function=_convertDensity

def _convertPressure(data):
    return data.convert("Pressure")
KnownRAMSESFields["Pressure"]._units=r"\rm{dyne}/\rm{cm}^{2}/\mu"
KnownRAMSESFields["Pressure"]._convert_function=_convertPressure

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = KnownRAMSESFields["%s-velocity" % ax]
    f._units = r"\rm{cm}/\rm{s}"
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
    if f not in KnownRAMSESFields:
        add_ramses_field(f, function=NullFunc, take_log=True,
                  validators = [ValidateDataField(f)],
                  particle_type = True)

for ax in 'xyz':
    KnownRAMSESFields["particle_velocity_%s" % ax]._convert_function = \
        _convertVelocity

def _convertParticleMass(data):
    return data.convert("mass")

KnownRAMSESFields["particle_mass"]._convert_function = \
        _convertParticleMass
KnownRAMSESFields["particle_mass"]._units = r"\mathrm{g}"

def _convertParticleMassMsun(data):
    return 1.0/mass_sun_cgs
add_field("ParticleMass", function=TranslationFunc("particle_mass"), 
          particle_type=True)
add_field("ParticleMassMsun",
          function=TranslationFunc("particle_mass"), 
          particle_type=True, convert_function=_convertParticleMassMsun)

def _Temperature(field, data):
    rv = data["Pressure"]/data["Density"]
    rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
    return rv
add_field("Temperature", function=_Temperature, units=r"\rm{K}")


# We now set up a couple particle fields.  This should eventually be abstracted
# into a single particle field function that adds them all on and is used
# across frontends, but that will need to wait until moving to using
# Coordinates, or vector fields.

def particle_count(field, data):
    pos = np.column_stack([data["particle_position_%s" % ax] for ax in 'xyz'])
    d = data.deposit(pos, method = "count")
    return d
RAMSESFieldInfo.add_field(("deposit", "%s_count" % "all"),
         function = particle_count,
         validators = [ValidateSpatial()],
         display_name = "\\mathrm{%s Count}" % "all",
         projection_conversion = '1')

def particle_mass(field, data):
    pos = np.column_stack([data["particle_position_%s" % ax] for ax in 'xyz'])
    d = data.deposit(pos, [data["ParticleMass"]], method = "sum")
    return d

RAMSESFieldInfo.add_field(("deposit", "%s_mass" % "all"),
         function = particle_mass,
         validators = [ValidateSpatial()],
         display_name = "\\mathrm{%s Mass}" % "all",
         units = r"\mathrm{g}",
         projected_units = r"\mathrm{g}\/\mathrm{cm}",
         projection_conversion = 'cm')

def particle_density(field, data):
    pos = np.column_stack([data["particle_position_%s" % ax] for ax in 'xyz'])
    d = data.deposit(pos, [data["ParticleMass"]], method = "sum")
    d /= data["CellVolume"]
    return d

RAMSESFieldInfo.add_field(("deposit", "%s_density" % "all"),
         function = particle_density,
         validators = [ValidateSpatial()],
         display_name = "\\mathrm{%s Density}" % "all",
         units = r"\mathrm{g}/\mathrm{cm}^{3}",
         projected_units = r"\mathrm{g}/\mathrm{cm}^{-2}",
         projection_conversion = 'cm')

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
             units = r"\rm{g}/\rm{cm}^3",
             projected_units = r"\rm{g}/\rm{cm}^2")
    add_field("%s_Fraction" % species,
             function=_SpeciesFraction,
             validators=ValidateDataField("%s_Density" % species),
             display_name="%s\/Fraction" % species)
    add_field("Comoving_%s_Density" % species,
             function=_SpeciesComovingDensity,
             validators=ValidateDataField("%s_Density" % species),
             display_name="Comoving\/%s\/Density" % species)
    add_field("%s_Mass" % species, units=r"\rm{g}", 
              function=_SpeciesMass, 
              validators=ValidateDataField("%s_Density" % species),
              display_name="%s\/Mass" % species)
    add_field("%s_MassMsun" % species, units=r"M_{\odot}", 
              function=_SpeciesMass, 
              convert_function=_convertCellMassMsun,
              validators=ValidateDataField("%s_Density" % species),
              display_name="%s\/Mass" % species)
    if _speciesMass.has_key(species):
        add_field("%s_NumberDensity" % species,
                  function=_SpeciesNumberDensity,
                  convert_function=_ConvertNumberDensity,
                  validators=ValidateDataField("%s_Density" % species))

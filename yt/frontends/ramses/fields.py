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
]

for f in known_ramses_particle_fields:
    if f not in KnownRAMSESFields:
        add_ramses_field(f, function=NullFunc, take_log=True,
                  validators = [ValidateDataField(f)],
                  particle_type = True)

def _ParticleMass(field, data):
    particles = data["particle_mass"].astype('float64') * \
                just_one(data["CellVolumeCode"].ravel())
    # Note that we mandate grid-type here, so this is okay
    return particles

def _convertParticleMass(data):
    return data.convert("Density")*(data.convert("cm")**3.0)
def _IOLevelParticleMass(grid):
    dd = dict(particle_mass = np.ones(1), CellVolumeCode=grid["CellVolumeCode"])
    cf = (_ParticleMass(None, dd) * _convertParticleMass(grid))[0]
    return cf
def _convertParticleMassMsun(data):
    return data.convert("Density")*((data.convert("cm")**3.0)/1.989e33)
def _IOLevelParticleMassMsun(grid):
    dd = dict(particle_mass = np.ones(1), CellVolumeCode=grid["CellVolumeCode"])
    cf = (_ParticleMass(None, dd) * _convertParticleMassMsun(grid))[0]
    return cf
add_field("ParticleMass",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMass,
          particle_convert_function=_IOLevelParticleMass)
add_field("ParticleMassMsun",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMassMsun,
          particle_convert_function=_IOLevelParticleMassMsun)


def _ParticleMass(field, data):
    particles = data["particle_mass"].astype('float64') * \
                just_one(data["CellVolumeCode"].ravel())
    # Note that we mandate grid-type here, so this is okay
    return particles

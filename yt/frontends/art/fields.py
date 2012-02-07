"""
ART-specific fields

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
    FieldInfo, \
    NullFunc, \
    TranslationFunc, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, mass_hydrogen_cgs

ARTFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = ARTFieldInfo.add_field

KnownARTFields = FieldInfoContainer()
add_art_field = KnownARTFields.add_field

translation_dict = {"Density":"density",
                    "TotalEnergy":"TotalEnergy",
                    "x-velocity":"velocity_x",
                    "y-velocity":"velocity_y",
                    "z-velocity":"velocity_z",
                    "Pressure":"pressure",
                    "Metallicity":"metallicity",
                    "GasEnergy":"GasEnergy",
                    "particle_mass":"ParticleMass"
                   }

for f,v in translation_dict.items():
    pfield = v.lower().startswith("particle")
    add_art_field(v, function=NullFunc, take_log=False,
                  validators = [ValidateDataField(v)],
                  particle_type = pfield)
    add_art_field(f, function=TranslationFunc(v), take_log=True,
                  particle_type = pfield)

#Particle Fields
def _get_convert(fname):
    def _conv(data):
        return 1.0
    return _conv

add_art_field("particle_mass", function=NullFunc, take_log=False,
              convert_function=_get_convert("particle_mass"),
              units=r"\rm{g}", particle_type=True)
    

def _convertDensity(data):
    return data.convert("Density")
KnownARTFields["Density"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTFields["Density"]._convert_function=_convertDensity

def _convertEnergy(data):
    return data.convert("GasEnergy")
KnownARTFields["GasEnergy"]._units = r"\rm{ergs}/\rm{g}"
KnownARTFields["GasEnergy"]._convert_function=_convertEnergy

def _Temperature(field, data):
    tr  = data["GasEnergy"] / data["Density"]
    tr /= data.pf.conversion_factors["GasEnergy"]
    tr *= data.pf.conversion_factors["Density"]
    return tr
def _convertTemperature(data):
    return data.convert("Temperature")
add_art_field("Temperature", function=_Temperature, units = r"\mathrm{K}")
KnownARTFields["Temperature"]._units = r"\mathrm{K}"
KnownARTFields["Temperature"]._convert_function=_convertTemperature

def _MetallicitySNII(field, data):
    #get the dimensionless mass fraction
    tr  = data["Metal_DensitySNII"] / data["Density"]
    tr *= data.pf.conversion_factors["Density"]    
    return tr
    
add_art_field("MetallicitySNII", function=_MetallicitySNII, units = r"\mathrm{K}")
KnownARTFields["MetallicitySNII"]._units = r"\mathrm{K}"

def _MetallicitySNIa(field, data):
    #get the dimensionless mass fraction
    tr  = data["Metal_DensitySNIa"] / data["Density"]
    tr *= data.pf.conversion_factors["Density"]    
    return tr
    
add_art_field("MetallicitySNIa", function=_MetallicitySNIa, units = r"\mathrm{K}")
KnownARTFields["MetallicitySNIa"]._units = r"\mathrm{K}"

def _Metallicity(field, data):
    #get the dimensionless mass fraction of the total metals
    tr  = data["Metal_DensitySNIa"] / data["Density"]
    tr += data["Metal_DensitySNII"] / data["Density"]
    tr *= data.pf.conversion_factors["Density"]    
    return tr
    
add_art_field("Metallicity", function=_Metallicity, units = r"\mathrm{K}")
KnownARTFields["Metallicity"]._units = r"\mathrm{K}"

def _Metal_Density(field,data):
    return data["Metal_DensitySNII"]+data["Metal_DensitySNIa"]
def _convert_Metal_Density(data):
    return data.convert("Metal_Density")

add_art_field("Metal_Density", function=_Metal_Density, units = r"\mathrm{K}")
KnownARTFields["Metal_Density"]._units = r"\mathrm{K}"
KnownARTFields["Metal_Density"]._convert_function=_convert_Metal_Density




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

KnownARTFields = FieldInfoContainer()
add_art_field = KnownARTFields.add_field

ARTFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = ARTFieldInfo.add_field

import numpy as na

#these are just the hydro fields
known_art_fields = [ 'Density','TotalEnergy',
                     'XMomentumDensity','YMomentumDensity','ZMomentumDensity',
                     'Pressure','Gamma','GasEnergy',
                     'MetalDensitySNII', 'MetalDensitySNIa',
                     'PotentialNew','PotentialOld']

#Add the fields, then later we'll individually defined units and names
for f in known_art_fields:
    if f not in ARTFieldInfo:
        add_field(f, function=lambda a,b: None, take_log=True,
                  validators = [ValidateDataField(f)])

#Hydro Fields that are verified to be OK unit-wise:
#Density
#Temperature

#Hydro Fields that need to be tested:
#TotalEnergy
#XYZMomentum
#Pressure
#Gamma
#GasEnergy
#MetalDensity SNII + SNia
#Potentials

#Hydro Derived fields that are untested:
#metallicities
#xyzvelocity

#Particle fields that are tested:
#particle_position_xyz
#particle_type
#particle_index
#particle_mass
#particle_mass_initial
#particle_age
#particle_velocity
#particle_metallicity12

#Particle fields that are untested:
#NONE


def _convertDensity(data):
    return data.convert("Density")
ARTFieldInfo["Density"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
ARTFieldInfo["Density"]._convert_function=_convertDensity

def _convertTotalEnergy(data):
    return data.convert("GasEnergy")
ARTFieldInfo["TotalEnergy"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["TotalEnergy"]._projected_units = r"\rm{K}"
ARTFieldInfo["TotalEnergy"]._convert_function=_convertTotalEnergy

def _convertXMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
ARTFieldInfo["XMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
ARTFieldInfo["XMomentumDensity"]._projected_units = r"\rm{K}"
ARTFieldInfo["XMomentumDensity"]._convert_function=_convertXMomentumDensity

def _convertYMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
ARTFieldInfo["YMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
ARTFieldInfo["YMomentumDensity"]._projected_units = r"\rm{K}"
ARTFieldInfo["YMomentumDensity"]._convert_function=_convertYMomentumDensity

def _convertZMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
ARTFieldInfo["ZMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
ARTFieldInfo["ZMomentumDensity"]._projected_units = r"\rm{K}"
ARTFieldInfo["ZMomentumDensity"]._convert_function=_convertZMomentumDensity

def _convertPressure(data):
    return data.convert("Pressure")
ARTFieldInfo["Pressure"]._units = r"\rm{g}/\rm{cm}/\rm{s}^2"
ARTFieldInfo["Pressure"]._projected_units = r"\rm{g}/\rm{s}^2"
ARTFieldInfo["Pressure"]._convert_function=_convertPressure

def _convertGamma(data):
    return 1.0
ARTFieldInfo["Gamma"]._units = r""
ARTFieldInfo["Gamma"]._projected_units = r""
ARTFieldInfo["Gamma"]._convert_function=_convertGamma

def _convertGasEnergy(data):
    return data.convert("GasEnergy")
ARTFieldInfo["GasEnergy"]._units = r"\rm{ergs}/\rm{g}"
ARTFieldInfo["GasEnergy"]._projected_units = r""
ARTFieldInfo["GasEnergy"]._convert_function=_convertGasEnergy

def _convertMetalDensitySNII(data):
    return data.convert("Density")
ARTFieldInfo["MetalDensitySNII"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["MetalDensitySNII"]._projected_units = r"\rm{g}/\rm{cm}^2"
ARTFieldInfo["MetalDensitySNII"]._convert_function=_convertMetalDensitySNII

def _convertMetalDensitySNIa(data):
    return data.convert("Density")
ARTFieldInfo["MetalDensitySNIa"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["MetalDensitySNIa"]._projected_units = r"\rm{g}/\rm{cm}^2"
ARTFieldInfo["MetalDensitySNIa"]._convert_function=_convertMetalDensitySNIa

def _convertPotentialNew(data):
    return data.convert("Potential")
ARTFieldInfo["PotentialNew"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["PotentialNew"]._projected_units = r"\rm{g}/\rm{cm}^2"
ARTFieldInfo["PotentialNew"]._convert_function=_convertPotentialNew

def _convertPotentialOld(data):
    return data.convert("Potential")
ARTFieldInfo["PotentialOld"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["PotentialOld"]._projected_units = r"\rm{g}/\rm{cm}^2"
ARTFieldInfo["PotentialOld"]._convert_function=_convertPotentialOld

####### Derived fields

def _temperature(field, data):
    tr  = data["GasEnergy"].astype('float64') #~1
    d = data["Density"].astype('float64')
    d[d==0.0] = -1.0 #replace the zeroes (that cause infs)
    tr /= d #
    assert na.all(na.isfinite(tr)) #diagnosing some problem...
    return tr
def _converttemperature(data):
    x  = data.pf.conversion_factors["Density"]
    x /= data.pf.conversion_factors["GasEnergy"]
    x *= data.pf.conversion_factors["Temperature"]
    return x
add_field("Temperature", function=_temperature, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Temperature"]._units = r"\mathrm{K}"
ARTFieldInfo["Temperature"]._projected_units = r"\mathrm{K}"
ARTFieldInfo["Temperature"]._convert_function=_converttemperature

def _metallicity_snII(field, data):
    tr  = data["MetalDensitySNII"] / data["Density"]
    return tr
add_field("Metallicity_SNII", function=_metallicity_snII, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Metallicity_SNII"]._units = r""
ARTFieldInfo["Metallicity_SNII"]._projected_units = r""

def _metallicity_snIa(field, data):
    tr  = data["MetalDensitySNIa"] / data["Density"]
    return tr
add_field("Metallicity_SNIa", function=_metallicity_snIa, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Metallicity_SNIa"]._units = r""
ARTFieldInfo["Metallicity_SNIa"]._projected_units = r""

def _x_velocity(data):
    tr  = data["XMomentumDensity"]/data["Density"]
    return tr
add_field("x_velocity", function=_x_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTFieldInfo["x_velocity"]._units = r"\rm{cm}/\rm{s}"
ARTFieldInfo["x_velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _y_velocity(data):
    tr  = data["YMomentumDensity"]/data["Density"]
    return tr
add_field("y_velocity", function=_y_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTFieldInfo["y_velocity"]._units = r"\rm{cm}/\rm{s}"
ARTFieldInfo["y_velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _z_velocity(data):
    tr  = data["ZMomentumDensity"]/data["Density"]
    return tr
add_field("z_velocity", function=_z_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTFieldInfo["z_velocity"]._units = r"\rm{cm}/\rm{s}"
ARTFieldInfo["z_velocity"]._projected_units = r"\rm{cm}/\rm{s}"


def _metal_density(field, data):
    tr  = data["MetalDensitySNIa"]
    tr += data["MetalDensitySNII"]
    return tr
add_field("Metal_Density", function=_metal_density, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Metal_Density"]._units = r""
ARTFieldInfo["Metal_Density"]._projected_units = r""


#Particle fields

#Derived particle fields


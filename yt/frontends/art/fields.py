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
import yt.utilities.lib as amr_utils
from yt.frontends.art.definitions import *

KnownARTFields = FieldInfoContainer()
add_art_field = KnownARTFields.add_field

ARTFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = ARTFieldInfo.add_field

import numpy as np

#these are just the hydro fields
#Add the fields, then later we'll individually defined units and names
for f in fluid_fields:
    add_art_field(f, function=NullFunc, take_log=True,
              validators = [ValidateDataField(f)])

for f in particle_fields:
    add_art_field(f, function=NullFunc, take_log=True,
              validators = [ValidateDataField(f)],
              particle_type = True)

add_art_field("particle_mass",function=NullFunc,take_log=True,
            validators=[ValidateDataField(f)],
            particle_type = True,
            convert_function= lambda x: x.convert("particle_mass"))

add_art_field("particle_mass_initial",function=NullFunc,take_log=True,
            validators=[ValidateDataField(f)],
            particle_type = True,
            convert_function= lambda x: x.convert("particle_mass"))

#Hydro Fields that are verified to be OK unit-wise:
#Density
#Temperature
#metallicities
#MetalDensity SNII + SNia

#Hydro Fields that need to be tested:
#TotalEnergy
#XYZMomentum
#Pressure
#Gamma
#GasEnergy
#Potentials
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

#Other checks:
#CellMassMsun == Density * CellVolume

def _convertDensity(data):
    return data.convert("Density")
KnownARTFields["Density"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTFields["Density"]._convert_function=_convertDensity

def _convertTotalEnergy(data):
    return data.convert("GasEnergy")
KnownARTFields["TotalEnergy"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["TotalEnergy"]._projected_units = r"\rm{K}"
KnownARTFields["TotalEnergy"]._convert_function=_convertTotalEnergy

def _convertXMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
KnownARTFields["XMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
KnownARTFields["XMomentumDensity"]._projected_units = r"\rm{K}"
KnownARTFields["XMomentumDensity"]._convert_function=_convertXMomentumDensity

def _convertYMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
KnownARTFields["YMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
KnownARTFields["YMomentumDensity"]._projected_units = r"\rm{K}"
KnownARTFields["YMomentumDensity"]._convert_function=_convertYMomentumDensity

def _convertZMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
KnownARTFields["ZMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
KnownARTFields["ZMomentumDensity"]._projected_units = r"\rm{K}"
KnownARTFields["ZMomentumDensity"]._convert_function=_convertZMomentumDensity

def _convertPressure(data):
    return data.convert("Pressure")
KnownARTFields["Pressure"]._units = r"\rm{g}/\rm{cm}/\rm{s}^2"
KnownARTFields["Pressure"]._projected_units = r"\rm{g}/\rm{s}^2"
KnownARTFields["Pressure"]._convert_function=_convertPressure

def _convertGamma(data):
    return 1.0
KnownARTFields["Gamma"]._units = r""
KnownARTFields["Gamma"]._projected_units = r""
KnownARTFields["Gamma"]._convert_function=_convertGamma

def _convertGasEnergy(data):
    return data.convert("GasEnergy")
KnownARTFields["GasEnergy"]._units = r"\rm{ergs}/\rm{g}"
KnownARTFields["GasEnergy"]._projected_units = r""
KnownARTFields["GasEnergy"]._convert_function=_convertGasEnergy

def _convertMetalDensitySNII(data):
    return data.convert('Density')
KnownARTFields["MetalDensitySNII"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["MetalDensitySNII"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTFields["MetalDensitySNII"]._convert_function=_convertMetalDensitySNII

def _convertMetalDensitySNIa(data):
    return data.convert('Density')
KnownARTFields["MetalDensitySNIa"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["MetalDensitySNIa"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTFields["MetalDensitySNIa"]._convert_function=_convertMetalDensitySNIa

def _convertPotentialNew(data):
    return data.convert("Potential")
KnownARTFields["PotentialNew"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["PotentialNew"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTFields["PotentialNew"]._convert_function=_convertPotentialNew

def _convertPotentialOld(data):
    return data.convert("Potential")
KnownARTFields["PotentialOld"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTFields["PotentialOld"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTFields["PotentialOld"]._convert_function=_convertPotentialOld

####### Derived fields

def _temperature(field, data):
    dg = data["GasEnergy"] #.astype('float64')
    dg /= data.pf.conversion_factors["GasEnergy"]
    dd = data["Density"] #.astype('float64')
    dd /= data.pf.conversion_factors["Density"]
    tr = dg/dd*data.pf.conversion_factors['tr']
    #ghost cells have zero density?
    tr[np.isnan(tr)] = 0.0
    #dd[di] = -1.0
    #if data.id==460:
    #tr[di] = -1.0 #replace the zero-density points with zero temp
    #print tr.min()
    #assert np.all(np.isfinite(tr))
    return tr
def _converttemperature(data):
    #x = data.pf.conversion_factors["Temperature"]
    x = 1.0
    return x
add_field("Temperature", function=_temperature, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Temperature"]._units = r"\mathrm{K}"
ARTFieldInfo["Temperature"]._projected_units = r"\mathrm{K}"
#ARTFieldInfo["Temperature"]._convert_function=_converttemperature

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

def _metallicity(field, data):
    tr  = data["Metal_Density"] / data["Density"]
    return tr
add_field("Metallicity", function=_metallicity, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Metallicity"]._units = r""
ARTFieldInfo["Metallicity"]._projected_units = r""

def _x_velocity(field,data):
    tr  = data["XMomentumDensity"]/data["Density"]
    return tr
add_field("x-velocity", function=_x_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTFieldInfo["x-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTFieldInfo["x-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _y_velocity(field,data):
    tr  = data["YMomentumDensity"]/data["Density"]
    return tr
add_field("y-velocity", function=_y_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTFieldInfo["y-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTFieldInfo["y-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _z_velocity(field,data):
    tr  = data["ZMomentumDensity"]/data["Density"]
    return tr
add_field("z-velocity", function=_z_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTFieldInfo["z-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTFieldInfo["z-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _metal_density(field, data):
    tr  = data["MetalDensitySNIa"]
    tr += data["MetalDensitySNII"]
    return tr
add_field("Metal_Density", function=_metal_density, units = r"\mathrm{K}",take_log=True)
ARTFieldInfo["Metal_Density"]._units = r""
ARTFieldInfo["Metal_Density"]._projected_units = r""


#Particle fields

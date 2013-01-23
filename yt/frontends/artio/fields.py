"""
ARTIO-specific fields

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
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
import numpy as np

KnownARTIOFields = FieldInfoContainer()
add_artio_field = KnownARTIOFields.add_field

#snl: doug removed RFI, but field name is needed in yt/data_objects/field_info_container.py?
ARTIOFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo, "RFI") 
add_field = ARTIOFieldInfo.add_field

known_artio_fields = [ 'Density','TotalEnergy',
                     'XMomentumDensity','YMomentumDensity','ZMomentumDensity',
                     'Pressure','Gamma','GasEnergy',
                     'MetalDensitySNII', 'MetalDensitySNIa',
                     'PotentialNew','PotentialOld']
                     
#Add the fields, then later we'll individually defined units and names
for f in known_artio_fields:
    add_artio_field(f, function=NullFunc, take_log=True,
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
KnownARTIOFields["Density"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["Density"]._convert_function=_convertDensity

def _convertTotalEnergy(data):
    return data.convert("GasEnergy")
KnownARTIOFields["TotalEnergy"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["TotalEnergy"]._projected_units = r"\rm{K}"
KnownARTIOFields["TotalEnergy"]._convert_function=_convertTotalEnergy

def _convertXMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
KnownARTIOFields["XMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
KnownARTIOFields["XMomentumDensity"]._projected_units = r"\rm{K}"
KnownARTIOFields["XMomentumDensity"]._convert_function=_convertXMomentumDensity

def _convertYMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
KnownARTIOFields["YMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
KnownARTIOFields["YMomentumDensity"]._projected_units = r"\rm{K}"
KnownARTIOFields["YMomentumDensity"]._convert_function=_convertYMomentumDensity

def _convertZMomentumDensity(data):
    tr  = data.convert("Mass")*data.convert("Velocity")
    tr *= (data.convert("Density")/data.convert("Mass"))
    return tr
KnownARTIOFields["ZMomentumDensity"]._units = r"\rm{mg}/\rm{s}/\rm{cm}^3"
KnownARTIOFields["ZMomentumDensity"]._projected_units = r"\rm{K}"
KnownARTIOFields["ZMomentumDensity"]._convert_function=_convertZMomentumDensity

def _convertPressure(data):
    return data.convert("Pressure")
KnownARTIOFields["Pressure"]._units = r"\rm{g}/\rm{cm}/\rm{s}^2"
KnownARTIOFields["Pressure"]._projected_units = r"\rm{g}/\rm{s}^2"
KnownARTIOFields["Pressure"]._convert_function=_convertPressure

def _convertGamma(data):
    return 1.0
KnownARTIOFields["Gamma"]._units = r""
KnownARTIOFields["Gamma"]._projected_units = r""
KnownARTIOFields["Gamma"]._convert_function=_convertGamma

def _convertGasEnergy(data):
    return data.convert("GasEnergy")
KnownARTIOFields["GasEnergy"]._units = r"\rm{ergs}/\rm{g}"
KnownARTIOFields["GasEnergy"]._projected_units = r""
KnownARTIOFields["GasEnergy"]._convert_function=_convertGasEnergy

def _convertMetalDensitySNII(data):
    return data.convert('Density')
KnownARTIOFields["MetalDensitySNII"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["MetalDensitySNII"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["MetalDensitySNII"]._convert_function=_convertMetalDensitySNII

def _convertMetalDensitySNIa(data):
    return data.convert('Density')
KnownARTIOFields["MetalDensitySNIa"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["MetalDensitySNIa"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["MetalDensitySNIa"]._convert_function=_convertMetalDensitySNIa

def _convertPotentialNew(data):
    return data.convert("Potential")
KnownARTIOFields["PotentialNew"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["PotentialNew"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["PotentialNew"]._convert_function=_convertPotentialNew

def _convertPotentialOld(data):
    return data.convert("Potential")
KnownARTIOFields["PotentialOld"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["PotentialOld"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["PotentialOld"]._convert_function=_convertPotentialOld

####### Derived fields

def _temperature(field, data):
    cd = data.pf.conversion_factors["Density"]
    cg = data.pf.conversion_factors["GasEnergy"]
    ct = data.pf.tr
    dg = data["GasEnergy"].astype('float64')
    dd = data["Density"].astype('float64')
    di = dd==0.0
    #dd[di] = -1.0
    tr = dg/dd
    #tr[np.isnan(tr)] = 0.0
    #if data.id==460:
    #    import pdb;pdb.set_trace()
    tr /= data.pf.conversion_factors["GasEnergy"]
    tr *= data.pf.conversion_factors["Density"]
    tr *= data.pf.tr
    #tr[di] = -1.0 #replace the zero-density points with zero temp
    #print tr.min()
    #assert np.all(np.isfinite(tr))
    return tr
def _converttemperature(data):
    x = data.pf.conversion_factors["Temperature"]
    x = 1.0
    return x
add_field("Temperature", function=_temperature, units = r"\mathrm{K}",take_log=True)
ARTIOFieldInfo["Temperature"]._units = r"\mathrm{K}"
ARTIOFieldInfo["Temperature"]._projected_units = r"\mathrm{K}"
ARTIOFieldInfo["Temperature"]._convert_function=_converttemperature

def _metallicity_snII(field, data):
    tr  = data["MetalDensitySNII"] / data["Density"]
    return tr
add_field("Metallicity_SNII", function=_metallicity_snII, units = r"\mathrm{K}",take_log=True)
ARTIOFieldInfo["Metallicity_SNII"]._units = r""
ARTIOFieldInfo["Metallicity_SNII"]._projected_units = r""

def _metallicity_snIa(field, data):
    tr  = data["MetalDensitySNIa"] / data["Density"]
    return tr
add_field("Metallicity_SNIa", function=_metallicity_snIa, units = r"\mathrm{K}",take_log=True)
ARTIOFieldInfo["Metallicity_SNIa"]._units = r""
ARTIOFieldInfo["Metallicity_SNIa"]._projected_units = r""

def _metallicity(field, data):
    tr  = data["Metal_Density"] / data["Density"]
    return tr
add_field("Metallicity", function=_metallicity, units = r"\mathrm{K}",take_log=True)
ARTIOFieldInfo["Metallicity"]._units = r""
ARTIOFieldInfo["Metallicity"]._projected_units = r""

def _x_velocity(data):
    tr  = data["XMomentumDensity"]/data["Density"]
    return tr
add_field("x-velocity", function=_x_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTIOFieldInfo["x-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTIOFieldInfo["x-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _y_velocity(data):
    tr  = data["YMomentumDensity"]/data["Density"]
    return tr
add_field("y-velocity", function=_y_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTIOFieldInfo["y-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTIOFieldInfo["y-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _z_velocity(data):
    tr  = data["ZMomentumDensity"]/data["Density"]
    return tr
add_field("z-velocity", function=_z_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTIOFieldInfo["z-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTIOFieldInfo["z-velocity"]._projected_units = r"\rm{cm}/\rm{s}"


def _metal_density(field, data):
    tr  = data["MetalDensitySNIa"]
    tr += data["MetalDensitySNII"]
    return tr
add_field("Metal_Density", function=_metal_density, units = r"\mathrm{K}",take_log=True)
ARTIOFieldInfo["Metal_Density"]._units = r""
ARTIOFieldInfo["Metal_Density"]._projected_units = r""


#Particle fields ######################

#Derived particle fields

def mass_dm(field, data):
    idx = data["particle_type"]<5
    #make a dumb assumption that the mass is evenly spread out in the grid
    #must return an array the shape of the grid cells
    tr  = data["Ones"] #create a grid in the right size
    if np.sum(idx)>0:
        tr /= np.prod(tr.shape) #divide by the volume
        tr *= np.sum(data['particle_mass'][idx]) #Multiply by total contaiend mass
        return tr
    else:
        return tr*0.0

add_field("particle_cell_mass_dm", function=mass_dm,
          validators=[ValidateSpatial(0)])


known_artio_particle_fields = [
    "POSITION_X",
    "POSITION_Y",
    "POSITION_Z",
    "VELOCITY_X",
    "VELOCITY_Y",
    "VELOCITY_Z",
    "particle_mass",
    "particle_mass_initial",
    "particle_creation_time"
    "particle_metallicitySNII"
    "particle_metallicitySNIA"
    "particle_identifier",
    "particle_refinement_level"
]
#    "particle_position_x",
#    "particle_position_y",
#    "particle_position_z",
#    "particle_velocity_x",
#    "particle_velocity_y",
#    "particle_velocity_z",

for f in known_artio_particle_fields:
    if f not in KnownARTIOFields:
        add_artio_field(f, function=NullFunc, take_log=True,
                  validators = [ValidateDataField(f)],
                  particle_type = True)

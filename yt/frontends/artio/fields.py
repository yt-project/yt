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
    TranslationFunc, \
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

known_artio_fields = ['Density', 'TotalEnergy',
                     'XMomentumDensity','YMomentumDensity','ZMomentumDensity',
                     'Pressure','GasEnergy',
                     'MetalDensitySNII', 'MetalDensitySNIa',
                     'Potential','PotentialHydro']
                     
#Add the fields, then later we'll individually defined units and names
for f in known_artio_fields:
    add_artio_field(f, function=NullFunc, take_log=True,
              validators = [ValidateDataField(f)])

add_artio_field("Gamma", function=NullFunc, take_log=False,
                validators = [ValidateDataField("Gamma")])

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

def _convertPotential(data):
    return data.convert("Potential")
KnownARTIOFields["Potential"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["Potential"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["Potential"]._convert_function=_convertPotential

def _convertPotentialHydro(data):
    return data.convert("Potential")
KnownARTIOFields["PotentialHydro"]._units = r"\rm{g}/\rm{cm}^3"
KnownARTIOFields["PotentialHydro"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownARTIOFields["PotentialHydro"]._convert_function=_convertPotentialHydro

####### Derived fields
import sys
def _temperature(field, data):
    tr = data["GasEnergy"]/data.pf.conversion_factors["GasEnergy"]*data.pf.conversion_factors["Density"]/data["Density"] #*(data["Gamma"]-1)*wmu 
    #ghost cells have zero density?
    tr[np.isnan(tr)] = 0.0
    #dd[di] = -1.0
    #if data.id==460:
    #tr[di] = -1.0 #replace the zero-density points with zero temp
    #print tr.min()
    #assert np.all(np.isfinite(tr))
    return tr
def _converttemperature(data):
    x = data.pf.conversion_factors["Temperature"]
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

def _x_velocity(field,data):
    tr  = data["XMomentumDensity"]/data["Density"]
    return tr
add_field("x-velocity", function=_x_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTIOFieldInfo["x-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTIOFieldInfo["x-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _y_velocity(field,data):
    tr  = data["YMomentumDensity"]/data["Density"]
    return tr
add_field("y-velocity", function=_y_velocity, units = r"\mathrm{cm/s}",take_log=False)
ARTIOFieldInfo["y-velocity"]._units = r"\rm{cm}/\rm{s}"
ARTIOFieldInfo["y-velocity"]._projected_units = r"\rm{cm}/\rm{s}"

def _z_velocity(field,data):
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

##################################################
#Particle fields

for ax in 'xyz':
    pf = "particle_velocity_%s" % ax
    add_artio_field(pf, function=NullFunc,
              particle_type=True)

add_artio_field("particle_mass", function=NullFunc, particle_type=True)

for ax in 'xyz':
    pf = "particle_position_%s" % ax
    add_artio_field(pf, function=NullFunc,
              particle_type=True)

def ParticleMass(field,data):
    return data['particle_mass']
add_field("ParticleMass",function=ParticleMass,units=r"\rm{g}",particle_type=True)


#Derived particle fields

def ParticleMassMsun(field,data):
    return data['particle_mass']/1.989e33 #*data.pf['Msun']
add_field("ParticleMassMsun",function=ParticleMassMsun,units=r"\rm{g}",particle_type=True)

def _creation_time(field,data):
    pa = data["particle_age"]
    tr = np.zeros(pa.shape,dtype='float')-1.0
    tr[pa>0] = pa[pa>0]
    return tr
add_field("creation_time",function=_creation_time,units=r"\rm{s}",particle_type=True)

def mass_dm(field, data):
    tr = np.ones(data.ActiveDimensions, dtype='float32')
    idx = data["particle_type"]<5
    #make a dumb assumption that the mass is evenly spread out in the grid
    #must return an array the shape of the grid cells
    if np.sum(idx)>0:
        tr /= np.prod(data['CellVolumeCode']*data.pf['mpchcm']**3.0) #divide by the volume
        tr *= np.sum(data['particle_mass'][idx])*data.pf['Msun'] #Multiply by total contaiend mass
        print tr.shape
        return tr
    else:
        return tr*1e-9

add_field("particle_cell_mass_dm", function=mass_dm, units = r"\mathrm{M_{sun}}",
        validators=[ValidateSpatial(0)],        
        take_log=False,
        projection_conversion="1")

def _spdensity(field, data):
    grid_mass = np.zeros(data.ActiveDimensions, dtype='float32')
    if data.star_mass.shape[0] ==0 : return grid_mass 
    amr_utils.CICDeposit_3(data.star_position_x,
                           data.star_position_y,
                           data.star_position_z,
                           data.star_mass.astype('float32'),
                           data.star_mass.shape[0],
                           grid_mass, 
                           np.array(data.LeftEdge).astype(np.float64),
                           np.array(data.ActiveDimensions).astype(np.int32), 
                           np.float64(data['dx']))
    return grid_mass 

#add_field("star_density", function=_spdensity,
#          validators=[ValidateSpatial(0)], convert_function=_convertDensity)

def _simple_density(field,data):
    mass = np.sum(data.star_mass)
    volume = data['dx']*data.ActiveDimensions.prod().astype('float64')
    return mass/volume

add_field("star_density", function=_simple_density,
          validators=[ValidateSpatial(0)], convert_function=_convertDensity)

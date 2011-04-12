"""
HaloProfiler re-centering functions.

Author: Stephen Skory <s@skory.us>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2011 Stephen Skory.  All Rights Reserved.

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

import numpy as na

from yt.funcs import *

from yt.data_objects.field_info_container import \
    add_field

centering_registry = {}

def add_function(name):
   def wrapper(func):
       centering_registry[name] = func
       return func
   return wrapper

#### Density ####

def _MatterDensityXTotalMass(field, data):
    return na.power((data['Dark_Matter_Density'] * data['TotalMassMsun']), 
                    1.)
def _Convert_MatterDensityXTotalMass(data):
    return 1
add_field("MatterDensityXTotalMass", units=r"",
          function=_MatterDensityXTotalMass,
          convert_function=_Convert_MatterDensityXTotalMass)

@add_function("Min_Dark_Matter_Density")
def find_minimum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('MatterDensityXTotalMass')
    return [mx,my,mz]

@add_function("Max_Dark_Matter_Density")
def find_maximum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('MatterDensityXTotalMass')
    return [mx,my,mz]

@add_function("CoM_Dark_Matter_Density")
def find_CoM_dm_density(data):
    dc_x = data.quantities['WeightedAverageQuantity']('x', 'MatterDensityXTotalMass')
    dc_y = data.quantities['WeightedAverageQuantity']('y', 'MatterDensityXTotalMass')
    dc_z = data.quantities['WeightedAverageQuantity']('z', 'MatterDensityXTotalMass')
    return [dc_x, dc_y, dc_z]

#### Temperature ####

def _TemperatureXTotalMass(field, data):
    return (data['Temperature'] * data['TotalMassMsun'])
def _Convert_TemperatureXTotalMass(data):
    return 1
add_field("TemperatureXTotalMass", units=r"",
          function=_TemperatureXTotalMass,
          convert_function=_Convert_TemperatureXTotalMass)

@add_function("Min_Temperature")
def find_minimum_temperature(data):
    ma, mini, mx, my, mz, mg = data.quantities['MaxLocation']('TemperatureXTotalMass')
    return [mx,my,mz]

@add_function("Max_Temperature")
def find_maximum_temperature(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('TemperatureXTotalMass')
    return [mx,my,mz]


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

#### Dark Matter Density ####

def _MatterDensityXDMMass(field, data):
    return (data['Dark_Matter_Density'] * data['Dark_Matter_Density'] \
        * data["CellVolume"])
def _Convert_MatterDensityXDMMass(data):
    return 1
add_field("MatterDensityXDMMass", units=r"",
          function=_MatterDensityXDMMass,
          convert_function=_Convert_MatterDensityXDMMass)

@add_function("Min_Dark_Matter_Density")
def find_minimum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('MatterDensityXDMMass')
    return [mx,my,mz]

@add_function("Max_Dark_Matter_Density")
def find_maximum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('MatterDensityXDMMass')
    return [mx,my,mz]

@add_function("CoM_Dark_Matter_Density")
def find_CoM_dm_density(data):
    dc_x = data.quantities['WeightedAverageQuantity']('x', 'MatterDensityXDMMass')
    dc_y = data.quantities['WeightedAverageQuantity']('y', 'MatterDensityXDMMass')
    dc_z = data.quantities['WeightedAverageQuantity']('z', 'MatterDensityXDMMass')
    return [dc_x, dc_y, dc_z]

#### Gas Density ####

def _GasDensityXCellMass(field, data):
    return (data['Density'] * data['CellMassMsun'])
def _Convert_GasDensityXCellMass(data):
    return 1
add_field("GasDensityXCellMass", units=r"",
          function=_GasDensityXCellMass,
          convert_function=_Convert_GasDensityXCellMass)

@add_function("Min_Gas_Density")
def find_minimum_gas_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('GasDensityXCellMass')
    return [mx,my,mz]

@add_function("Max_Gas_Density")
def find_maximum_gas_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('GasDensityXCellMass')
    return [mx,my,mz]

@add_function("CoM_Gas_Density")
def find_CoM_gas_density(data):
    dc_x = data.quantities['WeightedAverageQuantity']('x', 'GasDensityXCellMass')
    dc_y = data.quantities['WeightedAverageQuantity']('y', 'GasDensityXCellMass')
    dc_z = data.quantities['WeightedAverageQuantity']('z', 'GasDensityXCellMass')
    return [dc_x, dc_y, dc_z]

#### Total Density ####

def _TotalDensityXTotalMass(field, data):
    return (data['Density'] + data['Dark_Matter_Density']) * \
        data['TotalMassMsun']
def _Convert_TotalDensityXTotalMass(data):
    return 1
add_field("TotalDensityXTotalMass", units=r"",
          function=_TotalDensityXTotalMass,
          convert_function=_Convert_TotalDensityXTotalMass)

@add_function("Min_Total_Density")
def find_minimum_total_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('TotalDensityXTotalMass')
    return [mx,my,mz]

@add_function("Max_Total_Density")
def find_maximum_total_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('TotalDensityXTotalMass')
    return [mx,my,mz]

@add_function("CoM_Total_Density")
def find_CoM_total_density(data):
    dc_x = data.quantities['WeightedAverageQuantity']('x', 'TotalDensityXTotalMass')
    dc_y = data.quantities['WeightedAverageQuantity']('y', 'TotalDensityXTotalMass')
    dc_z = data.quantities['WeightedAverageQuantity']('z', 'TotalDensityXTotalMass')
    return [dc_x, dc_y, dc_z]

#### Temperature ####

def _TemperatureXCellMass(field, data):
    return (data['Temperature'] * data['CellMassMsun'])
def _Convert_TemperatureXCellMass(data):
    return 1
add_field("TemperatureXCellMass", units=r"",
          function=_TemperatureXCellMass,
          convert_function=_Convert_TemperatureXCellMass)

@add_function("Min_Temperature")
def find_minimum_temperature(data):
    ma, mini, mx, my, mz, mg = data.quantities['MinLocation']('TemperatureXCellMass')
    return [mx,my,mz]

@add_function("Max_Temperature")
def find_maximum_temperature(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('TemperatureXCellMass')
    return [mx,my,mz]


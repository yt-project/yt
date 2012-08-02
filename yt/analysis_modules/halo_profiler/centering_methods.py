"""
HaloProfiler re-centering functions.

Author: Stephen Skory <s@skory.us>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt-project.org/
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

@add_function("Min_Dark_Matter_Density")
def find_minimum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('Dark_Matter_Density',
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Dark_Matter_Density")
def find_maximum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Dark_Matter_Density',
                                                              preload=False)
    return (mx, my, mz)

@add_function("CoM_Dark_Matter_Density")
def find_CoM_dm_density(data):
   dc_x, dc_y, dc_z = data.quantities['CenterOfMass'](use_cells=False, 
                                                      use_particles=True,
                                                      preload=False)
   return (dc_x, dc_y, dc_z)

#### Gas Density ####

@add_function("Min_Gas_Density")
def find_minimum_gas_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('Density',
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Gas_Density")
def find_maximum_gas_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Density',
                                                              preload=False)
    return (mx, my, mz)

@add_function("CoM_Gas_Density")
def find_CoM_gas_density(data):
   dc_x, dc_y, dc_z = data.quantities['CenterOfMass'](use_cells=True, 
                                                      use_particles=False,
                                                      preload=False)
   return (dc_x, dc_y, dc_z)

#### Total Density ####

@add_function("Min_Total_Density")
def find_minimum_total_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('Matter_Density',
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Total_Density")
def find_maximum_total_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Matter_Density',
                                                              preload=False)
    return (mx, my, mz)

@add_function("CoM_Total_Density")
def find_CoM_total_density(data):
   dc_x, dc_y, dc_z = data.quantities['CenterOfMass'](use_cells=True, 
                                                      use_particles=True,
                                                      preload=False)
   return (dc_x, dc_y, dc_z)

#### Temperature ####

@add_function("Min_Temperature")
def find_minimum_temperature(data):
    ma, mini, mx, my, mz, mg = data.quantities['MinLocation']('Temperature',
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Temperature")
def find_maximum_temperature(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Temperature',
                                                              preload=False)
    return (mx, my, mz)


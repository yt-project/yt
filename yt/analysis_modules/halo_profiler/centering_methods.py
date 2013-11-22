"""
HaloProfiler re-centering functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

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
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Dark_Matter_Density")
def find_maximum_dm_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Dark_Matter_Density',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("CoM_Dark_Matter_Density")
def find_CoM_dm_density(data):
   dc_x, dc_y, dc_z = data.quantities['CenterOfMass'](use_cells=False, 
                                                      use_particles=True,
                                                      lazy_reader=True,
                                                      preload=False)
   return (dc_x, dc_y, dc_z)

#### Gas Density ####

@add_function("Min_Gas_Density")
def find_minimum_gas_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('Density',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Gas_Density")
def find_maximum_gas_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Density',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("CoM_Gas_Density")
def find_CoM_gas_density(data):
   dc_x, dc_y, dc_z = data.quantities['CenterOfMass'](use_cells=True, 
                                                      use_particles=False,
                                                      lazy_reader=True,
                                                      preload=False)
   return (dc_x, dc_y, dc_z)

#### Total Density ####

@add_function("Min_Total_Density")
def find_minimum_total_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MinLocation']('Matter_Density',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Total_Density")
def find_maximum_total_density(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Matter_Density',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("CoM_Total_Density")
def find_CoM_total_density(data):
   dc_x, dc_y, dc_z = data.quantities['CenterOfMass'](use_cells=True, 
                                                      use_particles=True,
                                                      lazy_reader=True,
                                                      preload=False)
   return (dc_x, dc_y, dc_z)

#### Temperature ####

@add_function("Min_Temperature")
def find_minimum_temperature(data):
    ma, mini, mx, my, mz, mg = data.quantities['MinLocation']('Temperature',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)

@add_function("Max_Temperature")
def find_maximum_temperature(data):
    ma, maxi, mx, my, mz, mg = data.quantities['MaxLocation']('Temperature',
                                                              lazy_reader=True,
                                                              preload=False)
    return (mx, my, mz)


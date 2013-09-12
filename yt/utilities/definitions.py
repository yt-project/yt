"""
Various definitions for various other modules and routines



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .physical_constants import \
   mpc_per_mpc, kpc_per_mpc, pc_per_mpc, au_per_mpc, rsun_per_mpc, \
   miles_per_mpc, km_per_mpc, cm_per_mpc, sec_per_Gyr, sec_per_Myr, \
   sec_per_year, sec_per_day

# The number of levels we expect to have at most
MAXLEVEL=48

axis_labels = [('y','z'),('x','z'),('x','y')]
axis_names = {0: 'x', 1: 'y', 2: 'z', 4:''}
inv_axis_names = {'x':0,'y':1,'z':2,
                  'X':0,'Y':1,'Z':2}

vm_axis_names = {0:'x', 1:'y', 2:'z', 3:'dx', 4:'dy'}

# The appropriate axes for which way we are slicing
x_dict = [1,0,0]
y_dict = [2,2,1]

x_names = ['y','x','x']
y_names = ['z','z','y']

# How many of each thing are in an Mpc
mpc_conversion = {'Mpc'   : mpc_per_mpc,
                  'mpc'   : mpc_per_mpc,
                  'kpc'   : kpc_per_mpc,
                  'pc'    : pc_per_mpc,
                  'au'    : au_per_mpc,
                  'rsun'  : rsun_per_mpc,
                  'miles' : miles_per_mpc,
                  'km'    : km_per_mpc,
                  'cm'    : cm_per_mpc}

# Nicely formatted versions of common length units
formatted_length_unit_names = {'mpc'     : 'Mpc',
                               'au'      : 'AU',
                               'rsun'    : 'R_\odot'}

# How many seconds are in each thing
sec_conversion = {'Gyr'   : sec_per_Gyr,
                  'Myr'   : sec_per_Myr,
                  'years' : sec_per_year,
                  'days'  : sec_per_day}

axis_labels = [('y','z'),('x','z'),('x','y')]

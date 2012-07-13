"""
Various definitions for various other modules and routines

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/

@todo: Move into yt.Defs, along with enki.EnkiDefs
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

from .physical_constants import \
   mpc_per_mpc, kpc_per_mpc, pc_per_mpc, au_per_mpc, rsun_per_mpc, \
   miles_per_mpc, cm_per_mpc, sec_per_Gyr, sec_per_Myr, sec_per_year, \
   sec_per_day

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
mpc_conversion = {'mpc'   : mpc_per_mpc,
                  'kpc'   : kpc_per_mpc,
                  'pc'    : pc_per_mpc,
                  'au'    : au_per_mpc,
                  'rsun'  : rsun_per_mpc,
                  'miles' : miles_per_mpc,
                  'cm'    : cm_per_mpc}

# How many seconds are in each thig
sec_conversion = {'Gyr'   : sec_per_Gyr,
                  'Myr'   : sec_per_Myr,
                  'years' : sec_per_year,
                  'days'  : sec_per_day}

axis_labels = [('y','z'),('x','z'),('x','y')]

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
   mpc_per_mpc, mpc_per_kpc, mpc_per_pc, mpc_per_au, mpc_per_rsun, \
   mpc_per_miles, mpc_per_cm, sec_per_Gyr, sec_per_Myr, sec_per_year, \
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
                  'kpc'   : mpc_per_kpc,
                  'pc'    : mpc_per_pc,
                  'au'    : mpc_per_au,
                  'rsun'  : mpc_per_rsun,
                  'miles' : mpc_per_miles,
                  'cm'    : mpc_per_cm}

# How many seconds are in each thig
sec_conversion = {'Gyr'   : sec_per_Gyr,
                  'Myr'   : sec_per_Myr,
                  'years' : sec_per_year,
                  'days'  : sec_per_day}

axis_labels = [('y','z'),('x','z'),('x','y')]

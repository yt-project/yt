"""
Halo quantity object.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Britton Smith, Matthew Turk.  All Rights Reserved.

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

from .halo_callbacks import HaloCallback
from .operator_registry import quantity_registry

import numpy as np

class HaloQuantity(HaloCallback):
    def __init__(self, quantity, function, *args, **kwargs):
        HaloCallback.__init__(self, function, args, kwargs)
        self.quantity = quantity
        
    def __call__(self, halo_catalog, halo):
        halo.quantities[self.quantity] = self.function(halo_catalog, halo, 
                                                       *self.args, **self.kwargs)
        return True

def add_quantity(name, function):
    quantity_registry[name] = HaloQuantity(name, function)

def center_of_mass(halo_catalog, halo):
    return (halo['particle_mass'] * 
            np.array([halo['particle_position_x'],
                      halo['particle_position_y'],
                      halo['particle_position_z']])).sum(axis=1) / \
                               halo['particle_mass'].sum()

add_quantity('center_of_mass', center_of_mass)

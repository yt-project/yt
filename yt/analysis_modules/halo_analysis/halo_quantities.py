"""
Halo quantity object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from .halo_callbacks import HaloCallback
from .operator_registry import quantity_registry

def add_quantity(name, function):
    quantity_registry[name] = HaloQuantity(name, function)

class HaloQuantity(HaloCallback):
    def __init__(self, function, *args, **kwargs):
        HaloCallback.__init__(self, function, args, kwargs)
        
    def __call__(self, halo):
        return self.function(halo_catalog, halo, 
                             *self.args, **self.kwargs)

def center_of_mass(halo):
    if halo.particles is None:
        raise RuntimeError("Center of mass requires halo to have particle data.")
    return (halo.particles['particle_mass'] * 
            np.array([halo.particles['particle_position_x'],
                      halo.particles['particle_position_y'],
                      halo.particles['particle_position_z']])).sum(axis=1) / \
                               halo.particles['particle_mass'].sum()

add_quantity('center_of_mass', center_of_mass)

def bulk_velocity(halo):
    if halo.particles is None:
        raise RuntimeError("Center of mass requires halo to have particle data.")
    return (halo.particles['particle_mass'] * 
            np.array([halo.particles['particle_velocity_x'],
                      halo.particles['particle_velocity_y'],
                      halo.particles['particle_velocity_z']])).sum(axis=1) / \
                               halo.particles['particle_mass'].sum()

add_quantity('bulk_velocity', bulk_velocity)

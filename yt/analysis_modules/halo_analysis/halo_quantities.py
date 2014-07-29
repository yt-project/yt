"""
Halo quantity object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.operator_registry import \
     OperatorRegistry

from .halo_callbacks import HaloCallback

quantity_registry = OperatorRegistry()

def add_quantity(name, function):
    quantity_registry[name] = HaloQuantity(function)

class HaloQuantity(HaloCallback):
    r"""
    A HaloQuantity is a function that takes minimally a Halo object, 
    performs some analysis, and then returns a value that is assigned 
    to an entry in the Halo.quantities dictionary.
    """
    def __init__(self, function, *args, **kwargs):
        HaloCallback.__init__(self, function, args, kwargs)
        
    def __call__(self, halo):
        return self.function(halo, *self.args, **self.kwargs)

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
        raise RuntimeError("Bulk velocity requires halo to have particle data.")
    return (halo.particles['particle_mass'] * 
            np.array([halo.particles['particle_velocity_x'],
                      halo.particles['particle_velocity_y'],
                      halo.particles['particle_velocity_z']])).sum(axis=1) / \
                               halo.particles['particle_mass'].sum()

add_quantity('bulk_velocity', bulk_velocity)

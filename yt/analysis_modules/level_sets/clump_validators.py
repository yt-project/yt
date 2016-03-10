"""
ClumpValidators and callbacks.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.lib.misc_utilities import \
    gravitational_binding_energy
from yt.utilities.operator_registry import \
    OperatorRegistry
from yt.utilities.physical_constants import \
    gravitational_constant_cgs as G

clump_validator_registry = OperatorRegistry()

def add_validator(name, function):
    clump_validator_registry[name] = ClumpValidator(function)

class ClumpValidator(object):
    r"""
    A ClumpValidator is a function that takes a clump and returns 
    True or False as to whether the clump is valid and shall be kept.
    """
    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None: self.args = []
        self.kwargs = kwargs
        if self.kwargs is None: self.kwargs = {}

    def __call__(self, clump):
        return self.function(clump, *self.args, **self.kwargs)
    
def _gravitationally_bound(clump, use_thermal_energy=True,
                           use_particles=True, truncate=True):
    "True if clump is gravitationally bound."

    use_particles &= \
      ("all", "particle_mass") in clump.data.ds.field_info
    
    bulk_velocity = clump.quantities.bulk_velocity(use_particles=use_particles)

    kinetic = 0.5 * (clump["gas", "cell_mass"] *
        ((bulk_velocity[0] - clump["gas", "velocity_x"])**2 +
         (bulk_velocity[1] - clump["gas", "velocity_y"])**2 +
         (bulk_velocity[2] - clump["gas", "velocity_z"])**2)).sum()

    if use_thermal_energy:
        kinetic += (clump["gas", "cell_mass"] *
                    clump["gas", "thermal_energy"]).sum()

    if use_particles:
        kinetic += 0.5 * (clump["all", "particle_mass"] *
            ((bulk_velocity[0] - clump["all", "particle_velocity_x"])**2 +
             (bulk_velocity[1] - clump["all", "particle_velocity_y"])**2 +
             (bulk_velocity[2] - clump["all", "particle_velocity_z"])**2)).sum()

    potential = clump.data.ds.quan(G *
        gravitational_binding_energy(
            clump["gas", "cell_mass"].in_cgs(),
            clump["index", "x"].in_cgs(),
            clump["index", "y"].in_cgs(),
            clump["index", "z"].in_cgs(),
            truncate, (kinetic / G).in_cgs()),
        kinetic.in_cgs().units)
    
    if truncate and potential >= kinetic:
        return True

    if use_particles:
        potential += clump.data.ds.quan(G *
            gravitational_binding_energy(
                clump["all", "particle_mass"].in_cgs(),
                clump["all", "particle_position_x"].in_cgs(),
                clump["all", "particle_position_y"].in_cgs(),
                clump["all", "particle_position_z"].in_cgs(),
                truncate, ((kinetic - potential) / G).in_cgs()),
        kinetic.in_cgs().units)

    return potential >= kinetic
add_validator("gravitationally_bound", _gravitationally_bound)

def _min_cells(clump, n_cells):
    "True if clump has a minimum number of cells."
    return (clump["index", "ones"].size >= n_cells)
add_validator("min_cells", _min_cells)

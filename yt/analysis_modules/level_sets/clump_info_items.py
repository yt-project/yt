"""
ClumpInfoCallback and callbacks.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.operator_registry import \
     OperatorRegistry

clump_info_registry = OperatorRegistry()

def add_clump_info(name, function):
    clump_info_registry[name] = ClumpInfoCallback(function)

class ClumpInfoCallback(object):
    r"""
    A ClumpInfoCallback is a function that takes a clump, computes a 
    quantity, and returns a string to be printed out for writing clump info.
    """
    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None: self.args = []
        self.kwargs = kwargs
        if self.kwargs is None: self.kwargs = {}

    def __call__(self, clump):
        return self.function(clump, *self.args, **self.kwargs)
    
def _total_cells(clump):
    n_cells = clump.data["index", "ones"].size
    return "Cells: %d." % n_cells
add_clump_info("total_cells", _total_cells)

def _cell_mass(clump):
    cell_mass = clump.data["gas", "cell_mass"].sum().in_units("Msun")
    return "Mass: %e Msun." % cell_mass
add_clump_info("cell_mass", _cell_mass)

def _mass_weighted_jeans_mass(clump):
    jeans_mass = clump.data.quantities.weighted_average_quantity(
        "jeans_mass", ("gas", "cell_mass")).in_units("Msun")
    return "Jeans Mass (mass-weighted): %.6e Msolar." % jeans_mass
add_clump_info("mass_weighted_jeans_mass", _mass_weighted_jeans_mass)

def _volume_weighted_jeans_mass(clump):
    jeans_mass = clump.data.quantities.weighted_average_quantity(
        "jeans_mass", ("index", "cell_volume")).in_units("Msun")
    return "Jeans Mass (volume-weighted): %.6e Msolar." % jeans_mass
add_clump_info("volume_weighted_jeans_mass", _volume_weighted_jeans_mass)

def _max_grid_level(clump):
    max_level = clump.data["index", "grid_level"].max()
    return "Max grid level: %d." % max_level
add_clump_info("max_grid_level", _max_grid_level)

def _min_number_density(clump):
    min_n = clump.data["gas", "number_density"].min().in_units("cm**-3")
    return "Min number density: %.6e cm^-3." % min_n
add_clump_info("min_number_density", _min_number_density)

def _max_number_density(clump):
    max_n = clump.data["gas", "number_density"].max().in_units("cm**-3")
    return "Max number density: %.6e cm^-3." % max_n
add_clump_info("max_number_density", _max_number_density)

def _distance_to_main_clump(clump, units="pc"):
    master = clump
    while master.parent is not None:
        master = master.parent
    master_com = clump.data.ds.arr(master.data.quantities.center_of_mass())
    my_com = clump.data.ds.arr(clump.data.quantities.center_of_mass())
    distance = np.sqrt(((master_com - my_com)**2).sum())
    return "Distance from master center of mass: %.6e %s." % \
      (distance.in_units(units), units)
add_clump_info("distance_to_main_clump", _distance_to_main_clump)

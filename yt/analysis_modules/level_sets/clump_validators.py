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

import numpy as np

from yt.utilities.operator_registry import \
     OperatorRegistry

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
    
def _gravitationally_bound(clump, truncate=True,
                           include_thermal_energy=True):
    "True if clump is gravitationally bound."
    return (clump.quantities.is_bound(truncate=truncate,
        include_thermal_energy=include_thermal_energy) > 1.0)
add_validator("gravitationally_bound", _gravitationally_bound)

def _min_cells(clump, n_cells):
    "True if clump has a minimum number of cells."
    return (clump["index", "ones"].size >= n_cells)
add_validator("min_cells", _min_cells)

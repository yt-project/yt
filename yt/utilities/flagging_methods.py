"""
Utilities for flagging zones for refinement in a dataset



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np # For modern purposes
from yt.utilities.lib.misc_utilities import grow_flagging_field
from yt.extern.six import add_metaclass

flagging_method_registry = {}

class RegisteredFlaggingMethod(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_type_name") and not cls._skip_add:
            flagging_method_registry[cls._type_name] = cls

@add_metaclass(RegisteredFlaggingMethod)
class FlaggingMethod(object):
    _skip_add = False

class OverDensity(FlaggingMethod):
    _type_name = "overdensity"
    def __init__(self, over_density):
        self.over_density = over_density

    def __call__(self, grid):
        rho = grid["density"] / (grid.ds.refine_by**grid.Level)
        return (rho > self.over_density)

class FlaggingGrid(object):
    def __init__(self, grid, methods):
        self.grid = grid
        flagged = np.zeros(grid.ActiveDimensions, dtype="bool")
        for method in methods:
            flagged |= method(self.grid)
        self.flagged = grow_flagging_field(flagged)
        self.subgrids = []
        self.left_index = grid.get_global_startindex()
        self.dimensions = grid.ActiveDimensions.copy()

    def find_subgrids(self):
        if not np.any(self.flagged): return []
        psg = ProtoSubgrid(self.flagged, self.left_index, self.dimensions)
        sgl = [psg]
        index = 0
        while index < len(sgl):
            psg = sgl[index]
            psg.shrink()
            if psg.dimensions.prod() == 0:
                sgl[index] = None
                continue
            while not psg.acceptable:
                new_psgs = []
                for i, dim in enumerate(np.argsort(psg.dimensions)[::-1]):
                    new_psgs = psg.find_by_zero_signature(dim)
                    if len(new_psgs) > 1:
                        break
                if len(new_psgs) <= 1:
                    new_psgs = psg.find_by_second_derivative()
                psg = new_psgs[0]
                sgl[index] = psg 
                sgl.extend(new_psgs[1:])
                psg.shrink()
            index += 1
        return sgl


# Much or most of this is directly translated from Enzo
class ProtoSubgrid(object):

    def __init__(self, flagged_base, left_index, dimensions, offset = (0,0,0)):
        self.left_index = left_index.copy()
        self.dimensions = dimensions.copy()
        self.flagged = flagged_base[offset[0]:offset[0]+dimensions[0],
                                    offset[1]:offset[1]+dimensions[1],
                                    offset[2]:offset[2]+dimensions[2]]
        self.compute_signatures()

    def compute_signatures(self):
        self.sigs = []
        for dim in range(3):
            d1 = (dim + 1) % 3
            d2 = int(dim == 0)
            self.sigs.append(self.flagged.sum(axis=d1).sum(axis=d2))

    @property
    def acceptable(self):
        return float(self.flagged.sum()) / self.flagged.size > 0.2

    def shrink(self):
        new_ind = []
        for dim in range(3):
            sig = self.sigs[dim]
            new_start = 0
            while sig[new_start] == 0:
                new_start += 1
            new_end = sig.size 
            while sig[new_end - 1] == 0:
                new_end -= 1
            self.dimensions[dim] = new_end - new_start
            self.left_index[dim] += new_start
            new_ind.append((new_start, new_end))
        self.flagged = self.flagged[new_ind[0][0]:new_ind[0][1],
                                    new_ind[1][0]:new_ind[1][1],
                                    new_ind[2][0]:new_ind[2][1]]
        self.compute_signatures()

    def find_by_zero_signature(self, dim):
        sig = self.sigs[dim]
        grid_ends = np.zeros((sig.size, 2), dtype='int64')
        ng = 0
        i = 0
        while i < sig.size:
            if sig[i] != 0:
                grid_ends[ng, 0] = i
                while i < sig.size and sig[i] != 0:
                    i += 1
                grid_ends[ng, 1] = i - 1
                ng += 1
            i += 1
        new_grids = []
        for si, ei in grid_ends[:ng,:]:
            li = self.left_index.copy()
            dims = self.dimensions.copy()
            li[dim] += si
            dims[dim] = ei - si
            offset = [0,0,0]
            offset[dim] = si
            new_grids.append(ProtoSubgrid(self.flagged, li, dims, offset))
        return new_grids

    def find_by_second_derivative(self):
        max_strength = 0
        max_axis = -1
        for dim in range(3):
            sig = self.sigs[dim]
            sd = sig[:-2] - 2.0*sig[1:-1] + sig[2:]
            center = int((self.flagged.shape[dim] - 1) / 2)
            strength = zero_strength = zero_cross = 0
            for i in range(1, sig.size-2):
                # Note that sd is offset by one
                if sd[i-1] * sd[i] < 0:
                    strength = np.abs(sd[i-1] - sd[i])
                    # TODO this differs from what I could find in ENZO
                    # there's |center - i| < |center - zero_cross| instead
                    # additionally zero_cross is undefined in first pass  
                    if strength > zero_strength or \
                       (strength == zero_strength and np.abs(center - i) < np.abs(zero_cross -i )):
                        zero_strength = strength
                        zero_cross = i
            if zero_strength > max_strength:
                max_axis = dim
        dims = self.dimensions.copy()
        li = self.left_index.copy()
        dims[max_axis] = zero_cross
        psg1 = ProtoSubgrid(self.flagged, li, dims)
        li[max_axis] += zero_cross
        dims[max_axis] = self.dimensions[max_axis] - zero_cross
        offset = np.zeros(3)
        offset[max_axis] = zero_cross
        psg2 = ProtoSubgrid(self.flagged, li, dims, offset)
        return [psg1, psg2]

    def __str__(self):
        return "LI: (%s) DIMS: (%s)" % (self.left_index, self.dimensions)

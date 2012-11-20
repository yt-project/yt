"""
Painting zones in a grid

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

import numpy as np

class FluidOperator(object):
    def apply(self, pf):
        for g in pf.h.grids: self(g)

class TopHatSphere(FluidOperator):
    def __init__(self, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        
    def __call__(self, grid, sub_select = None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax] - self.center[i])**2.0, r)
        np.sqrt(r, r)
        ind = (r <= self.radius)
        if sub_select is not None:
            ind &= sub_select
        for field, val in self.fields.iteritems():
            grid[field][ind] = val

class CoredSphere(FluidOperator):
    def __init__(self, core_radius, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        self.core_radius = core_radius

    def __call__(self, grid, sub_select = None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        r2 = self.radius**2
        cr2 = self.core_radius**2
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax] - self.center[i])**2.0, r)
        np.maximum(r, cr2, r)
        ind = (r <= r2)
        if sub_select is not None:
            ind &= sub_select
        for field, (outer_val, inner_val) in self.fields.iteritems():
            val = ((r[ind] - cr2) / (r2 - cr2))**0.5 * (outer_val - inner_val)
            grid[field][ind] = val + inner_val

class BetaModelSphere(FluidOperator):
    def __init__(self, beta, core_radius, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        self.core_radius = core_radius
        self.beta = beta
        
    def __call__(self, grid, sub_select = None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        r2 = self.radius**2
        cr2 = self.core_radius**2
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax] - self.center[i])**2.0, r)            
        ind = (r <= r2)
        if sub_select is not None:
            ind &= sub_select            
        for field, core_val in self.fields.iteritems() :
            val = core_val*(1.+r[ind]/cr2)**(-1.5*self.beta)
            grid[field][ind] = val

class NFWModelSphere(FluidOperator):
    def __init__(self, scale_radius, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        self.scale_radius = scale_radius
        
    def __call__(self, grid, sub_select = None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax] - self.center[i])**2.0, r)
        np,sqrt(r,r)
        ind = (r <= self.radius)
        r /= scale.radius
        if sub_select is not None:
            ind &= sub_select
        for field, scale_val in self.fields.iteritems() :
            val = scale_val/(r[ind]*(1.+r[ind])**2)
            grid[field][ind] = val
            
class RandomFluctuation(FluidOperator):
    def __init__(self, fields):
        self.fields = fields

    def __call__(self, grid, sub_select = None):
        if sub_select is None:
            sub_select = Ellipsis
        for field, mag in self.fields.iteritems():
            vals = grid[field][sub_select]
            rc = 1.0 + (np.random.random(vals.shape) - 0.5) * mag
            grid[field][sub_select] *= rc

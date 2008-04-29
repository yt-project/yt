"""
Quantities that can be derived from Enzo data that may also required additional
arguments.  (Standard arguments -- such as the center of a distribution of
points -- are excluded here, and left to the EnzoDerivedFields.)

Is this needlessly complicated?

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *

quantity_info = {}

class GridChildMaskWrapper:
    def __init__(self, grid):
        self.grid = grid
    def __getattr__(self, attr):
        return getattr(self.grid, attr)
    def __getitem__(self, item):
        return self.grid[item]*self.grid.child_mask

def func_wrapper(quantities_collection, quantities_object):
    func = quantities_object.function
    c_func = quantities_object.combine_function
    def call_func_unlazy(args, kwargs):
        retval = func(quantities_collection.data_object, *args, **kwargs)
        return c_func(quantities_collection.data_object, *retval)
    def call_func_lazy(args, kwargs):
        n_ret = quantities_object.n_ret
        retvals = [ [] for i in range(n_ret)]
        for g in quantities_collection.data_object._grids:
            rv = func(GridChildMaskWrapper(g), *args, **kwargs)
            for i in range(n_ret): retvals[i].append(rv[i])
            g.clear_data()
        retvals = [na.array(retvals[i]) for i in range(n_ret)]
        return c_func(quantities_collection.data_object, *retvals)
    def call_func(*args, **kwargs):
        lazy_reader = kwargs.pop('lazy_reader', False)
        if lazy_reader:
            return call_func_lazy(args, kwargs)
        else:
            return call_func_unlazy(args, kwargs)
    return call_func

class DerivedQuantity(object):
    def __init__(self, name, function,
                 combine_function, units = "",
                 n_ret = 0):
        self.name = name
        self.function = function
        self.combine_function = combine_function
        self.n_ret = n_ret

def add_quantity(name, **kwargs):
    if 'function' not in kwargs or 'combine_function' not in kwargs:
        mylog.error("Not adding field %s because both function and combine_function must be provided" % name)
        return
    f = kwargs.pop('function')
    c = kwargs.pop('combine_function')
    quantity_info[name] = DerivedQuantity(name, f, c, **kwargs)

class DerivedQuantityCollection(object):
    functions = quantity_info
    def __init__(self, data_object):
        self.data_object = data_object

    def __getitem__(self, key):
        if key not in self.functions:
            raise KeyError(key)
        return func_wrapper(self, self.functions[key])

def _CenterOfMass(data):
    x = (data["x"] * data["CellMassMsun"]).sum()
    y = (data["y"] * data["CellMassMsun"]).sum()
    z = (data["z"] * data["CellMassMsun"]).sum()
    den = data["CellMassMsun"].sum()
    return x,y,z, den
def _combCenterOfMass(data, x,y,z, den):
    return na.array([x.sum(), y.sum(), z.sum()])/den.sum()
add_quantity("CenterOfMass", function=_CenterOfMass,
             combine_function=_combCenterOfMass, n_ret = 4)

def _WeightedAverageQuantity(data, field, weight):
    num = (data[field] * data[weight]).sum()
    den = data[weight].sum()
    return num, den
def _combWeightedAverageQuantity(data, field, weight):
    return field.sum()/weight.sum()
add_quantity("WeightedAverageQuantity", function=_WeightedAverageQuantity,
             combine_function=_combWeightedAverageQuantity, n_ret = 2)

"""
# SpecificAngularMomentum is in cm^2 / s
j_vec = na.sum(sp["SpecificAngularMomentum"]*sp["CellMassMsun"],axis=1)/m_weight
j_mag = na.sqrt(na.sum(j_vec**2.0)) # cm^2 / s

eterm = na.sqrt(0.5*na.sum(sp["CellMassMsun"]*sp["VelocityMagnitude"]**2.0)/m_weight)

G = 6.67e-8 # cm^3 g^-1 s^-2
m_vir = na.sum(sp["CellMassMsun"]) + na.sum(sp["ParticleMassMsun"])# g
print "MVIR: %0.5e m_sun" % m_vir
m_vir *= 1.989e33
spin = j_mag * eterm / (m_vir * G)
"""

def _SpinParameter(data):
    weight = data["CellMassMsun"].sum()
    j_vec = (data["SpecificAngularMomentum"]*data["CellMassMsun"]).sum(axis=1)
    m_enc = data["CellMassMsun"].sum() + data["ParticleMassMsun"].sum()
    j_mag_pre = na.sum(j_vec**2.0)
    e_term_pre = 0.5*na.sum(data["CellMassMsun"]*data["VelocityMagnitude"]**2.0)
    return j_mag_pre, weight, m_enc, e_term_pre
def _combSpinParameter(data, j_mag_pre, weight, m_enc, e_term_pre):
    e_term = na.sqrt(e_term_pre.sum())/weight.sum()
    j_mag = na.sqrt(j_mag_pre.sum())/weight.sum()
    G = 6.67e-8 # cm^3 g^-1 s^-2
    spin = j_mag * e_term / (m_enc.sum() * G * 1.989e33)
    return spin
add_quantity("SpinParameter", function=_SpinParameter,
             combine_parameter=_combSpinParameter, n_ret=4)

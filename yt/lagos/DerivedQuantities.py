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
from yt.funcs import get_pbar

quantity_info = {}

class GridChildMaskWrapper:
    def __init__(self, grid, data_object):
        self.grid = grid
        self.data_object = data_object
    def __getattr__(self, attr):
        return getattr(self.grid, attr)
    def __getitem__(self, item):
        return self.data_object._get_data_from_grid(self.grid, item)

def func_wrapper(quantities_collection, quantities_object):
    func = quantities_object.function
    c_func = quantities_object.combine_function
    data_object = quantities_collection.data_object
    def call_func_unlazy(args, kwargs):
        retval = func(data_object, *args, **kwargs)
        return c_func(data_object, *retval)
    def call_func_lazy(args, kwargs):
        n_ret = quantities_object.n_ret
        retvals = [ [] for i in range(n_ret)]
        pbar = get_pbar("Calculating ", len(data_object._grids))
        for gi,g in enumerate(data_object._grids):
            rv = func(GridChildMaskWrapper(g, data_object), *args, **kwargs)
            for i in range(n_ret): retvals[i].append(rv[i])
            g.clear_data()
            pbar.update(gi)
        pbar.finish()
        retvals = [na.array(retvals[i]) for i in range(n_ret)]
        return c_func(data_object, *retvals)
    def call_func(*args, **kwargs):
        lazy_reader = kwargs.pop('lazy_reader', False)
        if lazy_reader and not quantities_object.force_unlazy:
            return call_func_lazy(args, kwargs)
        else:
            return call_func_unlazy(args, kwargs)
    return call_func

class DerivedQuantity(object):
    def __init__(self, name, function,
                 combine_function, units = "",
                 n_ret = 0, force_unlazy=False):
        self.name = name
        self.function = function
        self.combine_function = combine_function
        self.n_ret = n_ret
        self.force_unlazy = force_unlazy

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

    def keys(self):
        return self.functions.keys()

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

def _BulkVelocity(data):
    xv = (data["x-velocity"] * data["CellMassMsun"]).sum()
    yv = (data["y-velocity"] * data["CellMassMsun"]).sum()
    zv = (data["z-velocity"] * data["CellMassMsun"]).sum()
    w = data["CellMassMsun"].sum()
    return xv, yv, zv, w
def _combBulkVelocity(data, xv, yv, zv, w):
    w = w.sum()
    xv = xv.sum()/w
    yv = yv.sum()/w
    zv = zv.sum()/w
    return na.array([xv, yv, zv])
add_quantity("BulkVelocity", function=_BulkVelocity,
             combine_function=_combBulkVelocity, n_ret=4)

def _BaryonSpinParameter(data):
    am = data["SpecificAngularMomentum"]*data["CellMassMsun"]
    j_mag = am.sum(axis=1)
    m_enc = data["CellMassMsun"].sum() + data["ParticleMassMsun"].sum()
    e_term_pre = na.sum(data["CellMassMsun"]*data["VelocityMagnitude"]**2.0)
    weight=data["CellMassMsun"].sum()
    return j_mag, m_enc, e_term_pre, weight
def _combBaryonSpinParameter(data, j_mag, m_enc, e_term_pre, weight):
    # Because it's a vector field, we have to ensure we have enough dimensions
    if len(j_mag.shape) < 2: j_mag = na.expand_dims(j_mag, 0)
    W = weight.sum()
    M = m_enc.sum()
    J = na.sqrt(((j_mag.sum(axis=0))**2.0).sum())/W
    E = na.sqrt(e_term_pre.sum()/W)
    G = 6.67e-8 # cm^3 g^-1 s^-2
    spin = J * E / (M*1.989e33*G)
    print "WEIGHT", W, M, J, E, G, spin
    return spin
add_quantity("BaryonSpinParameter", function=_BaryonSpinParameter,
             combine_function=_combBaryonSpinParameter, n_ret=4)

def _IsBound(data, truncate = True):
    # Kinetic energy
    bv_x,bv_y,bv_z = data.quantities["BulkVelocity"]()
    kinetic = 0.5 * (data["CellMass"] * (
                       (data["x-velocity"] - bv_x)**2
                     + (data["y-velocity"] - bv_y)**2
                     + (data["z-velocity"] - bv_z)**2 )).sum()
    # Gravitational potential energy
    # We only divide once here because we have velocity in cgs, but radius is
    # in code.
    G = 6.67e-8 / data.convert("cm") # cm^3 g^-1 s^-2
    # the slow way (someone get clever!)
    total_cells = len(data['x'])
    cells_done = 0
    potential = 0.0
    pb = get_pbar("Calculating gravitational potential ", (0.5*(total_cells**2 - total_cells)))
    for q in xrange(total_cells-1):
        pot = data['CellMass'][(q+1):total_cells] / \
            (((data['x'][(q+1):total_cells]-data['x'][q])**2 + \
              (data['y'][(q+1):total_cells]-data['y'][q])**2 + \
              (data['z'][(q+1):total_cells]-data['z'][q])**2)**(0.5))
        potential += 2 * G * data['CellMass'][q] * pot.sum()
        cells_done += (total_cells - q - 1)
        pb.update(cells_done)
        if truncate and (potential > kinetic):
            pb.finish()
            return 1.0
    pb.finish()
    return (kinetic / potential)
def _combIsBound(data,bound):
    return bound
add_quantity("IsBound",function=_IsBound,combine_function=_combIsBound,n_ret=1,
             force_unlazy=True)

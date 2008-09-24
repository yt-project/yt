"""
Helpful recipes to make tasty AMR desserts

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

# Named imports
import yt.lagos as lagos
import yt.raven as raven
from yt.funcs import *
import numpy as na
import os.path, inspect, types
from functools import wraps
from yt.logger import ytLogger as mylog

def _fix_pf(pf):
    if isinstance(pf, lagos.StaticOutput): return pf
    if os.path.exists("%s.hierarchy" % pf):
        return lagos.EnzoStaticOutput(pf)
    elif pf.endswith(".hierarchy"):
        return lagos.EnzoStaticOutput(pf[:-10])
    # JS will have to implement the Orion one

__pf_centers = {}
def _fix_center(pf, center):
    if center is not None and iterable(center):
        center = na.array(center)
    else:
        if pf['CurrentTimeIdentifier'] in __pf_centers:
            center = __pf_centers[pf['CurrentTimeIdentifier']]
        else:
            center = pf.h.find_max("Density")[1]
            __pf_centers[pf['CurrentTimeIdentifier']] = center
    return center

def _fix_radius(pf, radius):
    if radius is not None:
        if iterable(radius):
            return radius[0] / pf[radius[1]]
        return radius
    mylog.info("Setting radius to be 0.1 of the box size")
    # yt-generalization : needs to be changed to 'unitary'
    return 0.1 / pf["1"]

def _fix_axis(pf, axis):
    if axis is None:
        raise ValueError("You need to specify an Axis!")
    elif isinstance(axis, types.IntType) and (0 <= x <= 2):
        return axis
    elif isinstance(axis, types.StringTypes) and \
         axis.upper() in 'XYZ':
        return 'XYZ'.find(axis.upper())
    else:
        raise ValueError("Invalid Axis specified.")


_arg_fixer = {
                'center' : (True, _fix_center),
                'radius' : (True, _fix_radius),
                'axis' : (True, _fix_axis),
             }

def fix_plot_args(plot_func):
    @wraps(plot_func)
    def call_func(self, pf, **kwargs):
        pf = _fix_pf(pf)
        fkwargs = inspect.getargspec(plot_func)[0]
        for arg in fkwargs:
            if arg in _arg_fixer:
                needs, fixer = _arg_fixer[arg]
                if arg in kwargs:
                    kwargs[arg] = fixer(pf, kwargs[arg])
                elif needs:
                    kwargs[arg] = fixer(pf, None)
        retval = plot_func(self, pf, **kwargs)
        if 'filename' in kwargs and 'filename' in fkwargs:
            retval.save(kwargs['filename'])
        return retval
    return call_func

which_pc = raven.PlotCollection

# to add:
#   zoom movies

class _RecipeBook(object):

    @fix_plot_args
    def mass_enclosed_radius(self, pf, center=None, radius=None, radius_unit="Radiuspc",
                             filename = None):
        pc = which_pc(pf, center=center)
        p = pc.add_profile_sphere(radius, '1', 
                                  [radius_unit, "CellMassMsun"],
                                  accumulation=True, weight=None)
        return pc

    @fix_plot_args
    def mass_enclosed_field(self, pf, field="Density", center=None, radius=None,
                            weight="CellMassMsun", filename = None):
        pc = which_pc(pf, center=center)
        p = pc.add_profile_sphere(radius, '1', 
                                  ["Radius", "CellMassMsun"],
                                  accumulation=True, weight=None)
        p.switch_z(field, weight=weight)
        p.switch_x("CellMassMsun")
        return pc

RecipeBook = _RecipeBook()

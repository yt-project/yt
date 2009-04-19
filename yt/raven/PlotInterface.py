"""
A high-level interface that obtains plots outside of a PlotCollection.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from yt.raven import *
from yt.funcs import *
from yt.recipes import _fix_pf

def check_plot_args(func):
    @wraps(func)
    def run_func(pf, *args, **kwargs):
        PlotTypes.Initialize()
        pf1 = _fix_pf(pf)
        if "center" not in kwargs:
            kwargs['center'] = pf1.h.find_max("Density")[1]
        p = func(pf1, *args, **kwargs)
        # This next bit is to ensure we have a strong reference
        # until the plot object is deleted
        if pf1 is not pf: p["pf"] = pf1
        return p
    return run_func

@check_plot_args
def get_slice(pf, field, axis, center=None, axes=None, figure=None,
              use_colorbar=True, **kwargs):
    """
    Get a single slice plot, with standard *field*, *axis* and *center*
    arguments.
    """
    coord = center[axis]
    data_source = pf.hierarchy.slice(axis, coord, field, center, **kwargs)
    plot = PlotTypes.SlicePlot(data_source, field,
            use_colorbar=use_colorbar, axes=axes, figure=figure)
    plot["Axis"] = lagos.axis_names[axis]
    return plot

@check_plot_args
def get_projection(pf, field, axis, weight_field=None, center=None,
                   axes=None, figure=None,
              use_colorbar=True, **kwargs):
    """
    Get a single projection plot, with standard *field*, *axis* and *center*
    arguments.
    """
    data_source = pf.hierarchy.proj(axis, field, weight_field, center=center, **kwargs)
    plot = PlotTypes.ProjectionPlot(data_source, field,
            use_colorbar=use_colorbar, axes=axes, figure=figure)
    plot["Axis"] = lagos.axis_names[axis]
    return plot



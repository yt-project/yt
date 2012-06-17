"""
Some convenience functions, objects, and iterators

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

import glob
import numpy as na
import os, os.path, inspect, types
from functools import wraps

# Named imports
from yt.funcs import *
from yt.config import ytcfg
from yt.utilities.parameter_file_storage import \
    output_type_registry, \
    simulation_time_series_registry, \
    EnzoRunDatabase

def load(*args ,**kwargs):
    """
    This function attempts to determine the base data type of a filename or
    other set of arguments by calling
    :meth:`yt.data_objects.api.StaticOutput._is_valid` until it finds a
    match, at which point it returns an instance of the appropriate
    :class:`yt.data_objects.api.StaticOutput` subclass.
    """
    if len(args) == 0:
        try:
            import Tkinter, tkFileDialog
        except ImportError:
            raise YTOutputNotIdentified(args, kwargs)
        root = Tkinter.Tk()
        filename = tkFileDialog.askopenfilename(parent=root,title='Choose a file')
        if filename != None:
            return load(filename)
        else:
            raise YTOutputNotIdentified(args, kwargs)
    candidates = []
    args = [os.path.expanduser(arg) if isinstance(arg, types.StringTypes)
            else arg for arg in args]
    valid_file = [os.path.exists(arg) if isinstance(arg, types.StringTypes) 
            else False for arg in args]
    if not any(valid_file):
        mylog.error("None of the arguments provided to load() is a valid file")
        mylog.error("Please check that you have used a correct path")
        raise YTOutputNotIdentified(args, kwargs)
    for n, c in output_type_registry.items():
        if n is None: continue
        if c._is_valid(*args, **kwargs): candidates.append(n)
    if len(candidates) == 1:
        return output_type_registry[candidates[0]](*args, **kwargs)
    if len(candidates) == 0:
        if ytcfg.get("yt", "enzo_db") != '' \
           and len(args) == 1 \
           and isinstance(args[0], types.StringTypes):
            erdb = EnzoRunDatabase()
            fn = erdb.find_uuid(args[0])
            n = "EnzoStaticOutput"
            if n in output_type_registry \
               and output_type_registry[n]._is_valid(fn):
                return output_type_registry[n](fn)
        mylog.error("Couldn't figure out output type for %s", args[0])
        raise YTOutputNotIdentified(args, kwargs)
    mylog.error("Multiple output type candidates for %s:", args[0])
    for c in candidates:
        mylog.error("    Possible: %s", c)
    raise YTOutputNotIdentified(args, kwargs)

def projload(pf, axis, weight_field = None):
    # This is something of a hack, so that we can just get back a projection
    # and not utilize any of the intermediate hierarchy objects.
    class ProjMock(dict):
        pass
    import h5py
    f = h5py.File(os.path.join(pf.fullpath, pf.parameter_filename + ".yt"))
    b = f["/Projections/%s/" % (axis)]
    wf = "weight_field_%s" % weight_field
    if wf not in b: raise KeyError(wf)
    fields = []
    for k in b:
        if k.startswith("weight_field"): continue
        if k.endswith("_%s" % weight_field):
            fields.append(k)
    proj = ProjMock()
    for f in ["px","py","pdx","pdy"]:
        proj[f] = b[f][:]
    for f in fields:
        new_name = f[:-(len(weight_field) + 1)]
        proj[new_name] = b[f][:]
    proj.axis = axis
    proj.pf = pf
    f.close()
    return proj

def simulation(parameter_filename, simulation_type):
    """
    Loads a simulation time series object of the specified
    simulation type.
    """

    if simulation_type not in simulation_time_series_registry:
        raise YTSimulationNotIdentified(simulation_type)

    return simulation_time_series_registry[simulation_type](parameter_filename)


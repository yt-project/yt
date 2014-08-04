"""
Some convenience functions, objects, and iterators



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import numpy as np
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
    :meth:`yt.data_objects.api.Dataset._is_valid` until it finds a
    match, at which point it returns an instance of the appropriate
    :class:`yt.data_objects.api.Dataset` subclass.
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
    valid_file = []
    for argno, arg in enumerate(args):
        if isinstance(arg, types.StringTypes):
            if os.path.exists(arg):
                valid_file.append(True)
            elif arg.startswith("http"):
                valid_file.append(True)
            else:
                if os.path.exists(os.path.join(ytcfg.get("yt", "test_data_dir"), arg)):
                    valid_file.append(True)
                    args[argno] = os.path.join(ytcfg.get("yt", "test_data_dir"), arg)
                else:
                    valid_file.append(False)
        else:
            valid_file.append(False)
    if not any(valid_file):
        try:
            from yt.data_objects.time_series import DatasetSeries
            ts = DatasetSeries.from_filenames(*args, **kwargs)
            return ts
        except YTOutputNotIdentified:
            pass
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
            n = "EnzoDataset"
            if n in output_type_registry \
               and output_type_registry[n]._is_valid(fn):
                return output_type_registry[n](fn)
        mylog.error("Couldn't figure out output type for %s", args[0])
        raise YTOutputNotIdentified(args, kwargs)
    mylog.error("Multiple output type candidates for %s:", args[0])
    for c in candidates:
        mylog.error("    Possible: %s", c)
    raise YTOutputNotIdentified(args, kwargs)

def projload(ds, axis, weight_field = None):
    # This is something of a hack, so that we can just get back a projection
    # and not utilize any of the intermediate index objects.
    class ProjMock(dict):
        pass
    import h5py
    f = h5py.File(os.path.join(ds.fullpath, ds.parameter_filename + ".yt"))
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
    proj.ds = ds
    f.close()
    return proj

def simulation(parameter_filename, simulation_type, find_outputs=False):
    """
    Loads a simulation time series object of the specified
    simulation type.
    """

    if simulation_type not in simulation_time_series_registry:
        raise YTSimulationNotIdentified(simulation_type)

    if os.path.exists(parameter_filename):
        valid_file = True
    elif os.path.exists(os.path.join(ytcfg.get("yt", "test_data_dir"),
                                     parameter_filename)):
        parameter_filename = os.path.join(ytcfg.get("yt", "test_data_dir"),
                                          parameter_filename)
        valid_file = True
    else:
        valid_file = False
        
    if not valid_file:
        raise YTOutputNotIdentified((parameter_filename, simulation_type), 
                                    dict(find_outputs=find_outputs))
    
    return simulation_time_series_registry[simulation_type](parameter_filename, 
                                                            find_outputs=find_outputs)


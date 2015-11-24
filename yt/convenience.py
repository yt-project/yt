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

import os

# Named imports
from yt.extern.six import string_types
from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.parameter_file_storage import \
    output_type_registry, \
    simulation_time_series_registry, \
    EnzoRunDatabase
from yt.utilities.exceptions import \
    YTOutputNotIdentified, \
    YTSimulationNotIdentified
from yt.utilities.hierarchy_inspection import find_lowest_subclasses

def load(*args ,**kwargs):
    """
    This function attempts to determine the base data type of a filename or
    other set of arguments by calling
    :meth:`yt.data_objects.api.Dataset._is_valid` until it finds a
    match, at which point it returns an instance of the appropriate
    :class:`yt.data_objects.api.Dataset` subclass.
    """
    candidates = []
    args = [os.path.expanduser(arg) if isinstance(arg, string_types)
            else arg for arg in args]
    valid_file = []
    for argno, arg in enumerate(args):
        if isinstance(arg, string_types):
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
    types_to_check = output_type_registry
    if not any(valid_file):
        try:
            from yt.data_objects.time_series import DatasetSeries
            ts = DatasetSeries.from_filenames(*args, **kwargs)
            return ts
        except (TypeError, YTOutputNotIdentified):
            pass
        # We check if either the first argument is a dict or list, in which
        # case we try identifying candidates.
        if len(args) > 0 and isinstance(args[0], (list, dict)):
            # This fixes issues where it is assumed the first argument is a
            # file
            types_to_check = dict((n, v) for n, v in
                    output_type_registry.items() if n.startswith("stream_"))
            # Better way to do this is to override the output_type_registry
        else:
            mylog.error("None of the arguments provided to load() is a valid file")
            mylog.error("Please check that you have used a correct path")
            raise YTOutputNotIdentified(args, kwargs)
    for n, c in types_to_check.items():
        if n is None: continue
        if c._is_valid(*args, **kwargs): candidates.append(n)

    # convert to classes
    candidates = [output_type_registry[c] for c in candidates]
    # Find only the lowest subclasses, i.e. most specialised front ends
    candidates = find_lowest_subclasses(candidates)
    if len(candidates) == 1:
        return candidates[0](*args, **kwargs)
    if len(candidates) == 0:
        if ytcfg.get("yt", "enzo_db") != '' \
           and len(args) == 1 \
           and isinstance(args[0], string_types):
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


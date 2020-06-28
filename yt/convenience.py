import os

# Named imports
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

def _sanitize_load_args(*args):
    """Filter out non-pathlike arguments, ensure list form, and expand '~' tokens"""
    try:
        # os.PathLike is python >= 3.6
        path_types = str, os.PathLike
    except AttributeError:
        path_types = str,

    return [os.path.expanduser(arg) if isinstance(arg, path_types)
            else arg for arg in args]

def load(*args ,**kwargs):
    """
    This function attempts to determine the base data type of a filename or
    other set of arguments by calling
    :meth:`yt.data_objects.static_output.Dataset._is_valid` until it finds a
    match, at which point it returns an instance of the appropriate
    :class:`yt.data_objects.static_output.Dataset` subclass.
    """
    args = _sanitize_load_args(*args)

    if any([wildcard in args[0] for wildcard in "[]?!*"]):
        from yt.data_objects.time_series import DatasetSeries
        return DatasetSeries(*args, **kwargs)

    if not (os.path.exists(args[0]) or args[0].startswith("http")):
        test_path = os.path.join(ytcfg.get("yt", "test_data_dir"), args[0])
        if os.path.exists(test_path):
            args[0] = test_path
        else:
            raise OSError("No such file or directory: %s" % args[0])

    types_to_check = output_type_registry
    candidates = []
    for n, c in types_to_check.items():
        if n is not None and c._is_valid(*args, **kwargs):
            candidates.append(n)

    # convert to classes
    candidates = [output_type_registry[c] for c in candidates]
    # Find only the lowest subclasses, i.e. most specialised front ends
    candidates = find_lowest_subclasses(candidates)
    if len(candidates) == 1:
        return candidates[0](*args, **kwargs)
    if len(candidates) == 0:
        if ytcfg.get("yt", "enzo_db") != '' \
           and len(args) == 1 \
           and isinstance(args[0], str):
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

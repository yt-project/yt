import os

# Named imports
from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.parameter_file_storage import \
    output_type_registry, \
    simulation_time_series_registry
from yt.utilities.exceptions import \
    YTOutputNotIdentified, \
    YTSimulationNotIdentified
from yt.utilities.hierarchy_inspection import find_lowest_subclasses

def load(fn, *args, **kwargs):
    """
    This function attempts to determine the base data type of a filename or
    other set of arguments by calling
    :meth:`yt.data_objects.static_output.Dataset._is_valid` until it finds a
    match, at which point it returns an instance of the appropriate
    :class:`yt.data_objects.static_output.Dataset` subclass.
    """
    fn = os.path.expanduser(fn)

    if any([wildcard in fn for wildcard in "[]?!*"]):
        from yt.data_objects.time_series import DatasetSeries
        return DatasetSeries(fn, *args, **kwargs)

    if not (os.path.exists(fn) or fn.startswith("http")):
        test_path = os.path.join(ytcfg.get("yt", "test_data_dir"), fn)
        if os.path.exists(test_path):
            fn = test_path
        else:
            raise OSError("No such file or directory: %s" % fn)

    candidates = []
    for cls in output_type_registry.values():
        if cls._is_valid(fn, *args, **kwargs):
            candidates.append(cls)

    # Find only the lowest subclasses, i.e. most specialised front ends
    candidates = find_lowest_subclasses(candidates)

    if len(candidates) == 1:
        return candidates[0](fn, *args, **kwargs)

    if len(candidates) > 1:
        mylog.error("Multiple output type candidates for %s:", fn)
        for c in candidates:
            mylog.error("    Possible: %s", c)

    raise YTOutputNotIdentified([fn, *args], kwargs)

def simulation(fn, simulation_type, find_outputs=False):
    """
    Loads a simulation time series object of the specified
    simulation type.
    """

    if not os.path.exists(fn):
        test_path = os.path.join(ytcfg.get("yt", "test_data_dir"), fn)
        if os.path.exists(test_path):
            fn = test_path
        else:
            raise OSError("No such file or directory: %s" % fn)

    try:
        cls = simulation_time_series_registry[simulation_type]
    except KeyError:
        raise YTSimulationNotIdentified(simulation_type)

    return cls(fn, find_outputs=find_outputs)

def load_enzo_db(fn):
    from yt.utilities.parameter_file_storage import EnzoRunDatabase
    from yt.frontends.enzo import EnzoDataset
    if not ytcfg.get("yt", "enzo_db"):
        raise OSError("enzo_db location is not properly setup.")

    erdb = EnzoRunDatabase()
    fn = os.path.expanduser(fn)
    fn = erdb.find_uuid(fn)

    if not EnzoDataset._is_valid(fn):
        raise YTOutputNotIdentified(fn, {})

    return EnzoDataset(fn)

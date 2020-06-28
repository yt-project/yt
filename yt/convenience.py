import os

# Named imports
from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.exceptions import YTOutputNotIdentified, YTSimulationNotIdentified
from yt.utilities.hierarchy_inspection import find_lowest_subclasses
from yt.utilities.parameter_file_storage import (
    output_type_registry,
    simulation_time_series_registry,
)


def load(fn, *args, **kwargs):
    """
    Load a Dataset or DatasetSeries object.
    The data format is automatically discovered, and the exact return type is the
    corresponding subclass of :class:`yt.data_objects.static_output.Dataset`.
    A :class:`yt.data_objects.time_series.DatasetSeries` is created if the first
    argument is a pattern.

    Parameters
    ----------
    fn : str, os.Pathlike, or byte (types supported by os.path.expandusers)
        A path to the data location. This can be a file name, directory name, a glob
        pattern, or a url (for data types that support it).

    Additional arguments, if any, are passed down to the return class.

    Returns
    -------
    :class:`yt.data_objects.static_output.Dataset` object
        If fn is a single path, create a Dataset from the appropriate subclass.

    :class:`yt.data_objects.time_series.DatasetSeries`
        If fn is a glob pattern (i.e. containing wildcards '[]?!*'), create a series.

    Raises
    ------
    OSError
        If fn does not match any existing file or directory.

    yt.utilities.exceptions.YTOutputNotIdentified
        If fn matches existing files or directories with undetermined format.
    """
    fn = os.path.expanduser(fn)

    if any(wildcard in fn for wildcard in "[]?!*"):
        from yt.data_objects.time_series import DatasetSeries

        return DatasetSeries(fn, *args, **kwargs)

    if not (os.path.exists(fn) or fn.startswith("http")):
        data_dir = ytcfg.get("yt", "test_data_dir")
        alt_fn = os.path.join(data_dir, fn)
        if os.path.exists(alt_fn):
            fn = alt_fn
        else:
            msg = "No such file or directory: %s" % fn
            if os.path.exists(data_dir):
                msg += "\n(Also tried %s)" % alt_fn
            raise OSError(msg)

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
    Load a simulation time series object of the specified simulation type.

    Parameters
    ----------
    fn : str, os.Pathlike, or byte (types supported by os.path.expandusers)
        Name of the data file or directory.

    simulation_type : str
        E.g. 'Enzo'

    find_outputs : bool
        Defaults to False

    Raises
    ------
    OSError
        If fn is not found.

    yt.utilities.exceptions.YTSimulationNotIdentified
        If simulation_type is unknown.
    """

    if not os.path.exists(fn):
        alt_fn = os.path.join(ytcfg.get("yt", "test_data_dir"), fn)
        if os.path.exists(alt_fn):
            fn = alt_fn
        else:
            raise OSError("No such file or directory: %s" % fn)

    try:
        cls = simulation_time_series_registry[simulation_type]
    except KeyError:
        raise YTSimulationNotIdentified(simulation_type)

    return cls(fn, find_outputs=find_outputs)

from yt.config import ytcfg
from yt.utilities.on_demand_imports import _dask as dask


def is_delayed(obj):
    """checks if obj is a dask array (or a unyt-dask array)"""
    return isinstance(obj, dask.array.Array)


def _empty_decorator(func):
    def call_func(*args, **kwargs):
        return func(*args, **kwargs)

    return call_func


def _passthrough_function(*args, **kwargs):
    return args


def _get_dask_delayed(use_dask: bool):
    if use_dask and ytcfg.get("yt", "internals", "dask_enabled"):
        return dask.delayed  # contains additional internal config checks
    else:
        return _empty_decorator


def _get_dask_compute(use_dask: bool):
    if use_dask and ytcfg.get("yt", "internals", "dask_enabled"):
        return dask.compute
    else:
        return _passthrough_function

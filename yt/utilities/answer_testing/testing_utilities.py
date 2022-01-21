import functools
import hashlib
import inspect
import os

import numpy as np
import pytest
import yaml

import yt.visualization.particle_plots as particle_plots
import yt.visualization.plot_window as pw
import yt.visualization.profile_plotter as profile_plotter
from yt.config import ytcfg
from yt.data_objects.selection_objects.region import YTRegion
from yt.data_objects.static_output import Dataset
from yt.loaders import load, load_simulation
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.visualization.volume_rendering.scene import Scene


def _streamline_for_io(params):
    r"""
    Put test results in a more io-friendly format.

    Many of yt's tests use objects such as tuples as test parameters
    (fields, for instance), but when these objects are written to a
    yaml file, yaml includes python specific anchors that make the file
    harder to read and less portable. The goal of this function is to
    convert these objects to strings (using __repr__() has it's own
    issues) in order to solve this problem.

    Parameters
    ----------
    params : dict
        The dictionary of test parameters in the form
        {param_name : param_value}.

    Returns
    -------
    streamlined_params : dict
        The dictionary of parsed and converted
        {param_name : param_value} pairs.
    """
    streamlined_params = {}
    for key, value in params.items():
        # Check for user-defined functions
        if inspect.isfunction(key):
            key = key.__name__
        if inspect.isfunction(value):
            value = value.__name__
        # The key can be nested iterables, e.g.,
        # d = [None, ('sphere', (center, (0.1, 'unitary')))] so we need
        # to use recursion
        if not isinstance(key, str) and hasattr(key, "__iter__"):
            key = _iterable_to_string(key)
        # The value can also be nested iterables
        if not isinstance(value, str) and hasattr(value, "__iter__"):
            value = _iterable_to_string(value)
        # Scene objects need special treatment to make them more IO friendly
        if isinstance(value, Scene):
            value = "Scene"
        elif isinstance(value, YTRegion):
            value = "Region"
        streamlined_params[key] = value
    return streamlined_params


def _iterable_to_string(iterable):
    r"""
    An extension of streamline_for_io that does the work of making an
    iterable more io-friendly.

    Parameters
    ----------
    iterable : python iterable
        The object to be parsed and converted.

    Returns
    -------
    result : str
        The io-friendly version of the passed iterable.
    """
    result = iterable.__class__.__name__
    for elem in iterable:
        # Check for user-defined functions
        if inspect.isfunction(elem):
            result += "_" + elem.__name__
        # Non-string iterables (e.g., lists, tuples, etc.)
        elif not isinstance(elem, str) and hasattr(elem, "__iter__"):
            result += "_" + _iterable_to_string(elem)
        # Non-string non-iterables (ints, floats, etc.)
        elif not isinstance(elem, str) and not hasattr(elem, "__iter__"):
            result += "_" + str(elem)
        # Strings
        elif isinstance(elem, str):
            result += "_" + elem
    return result


def _hash_results(results):
    r"""
    Driver function for hashing the test result.

    Parameters
    ----------
    results : dict
        Dictionary of {test_name : test_result} pairs.

    Returns
    -------
    results : dict
        Same as the passed results, but the test_results are now
        hex digests of the hashed test_result.
    """
    # Here, results should be comprised of only the tests, not the test
    # parameters
    # Use a new dictionary so as to not overwrite the non-hashed test
    # results in case those are to be saved
    hashed_results = {}
    for test_name, test_value in results.items():
        hashed_results[test_name] = generate_hash(test_value)
    return hashed_results


def _hash_dict(data):
    r"""
    Specifically handles hashing a dictionary object.

    Parameters
    ----------
    data : dict
        The dictionary to be hashed.

    Returns
    -------
    hd.hexdigest : str
        The hex digest of the hashed dictionary.
    """
    hd = None
    for key, value in data.items():
        # Some keys are tuples, not strings
        if not isinstance(key, str):
            key = key.__repr__()
        # Test suites can return values that are dictionaries of other tests
        if isinstance(value, dict):
            hashed_data = _hash_dict(value)
        else:
            hashed_data = bytearray(key.encode("utf8")) + bytearray(value)
        # If we're returning from a recursive call (and therefore hashed_data
        # is a hex digest), we need to encode the string before it can be
        # hashed
        if isinstance(hashed_data, str):
            hashed_data = hashed_data.encode("utf8")
        if hd is None:
            hd = hashlib.md5(hashed_data)
        else:
            hd.update(hashed_data)
    return hd.hexdigest()


def generate_hash(data):
    r"""
    Actually performs the hash operation.

    Parameters
    ----------
    data : python object
        Data to be hashed.

    Returns
    -------
    hd : str
        Hex digest of the hashed data.
    """
    # Sometimes md5 complains that the data is not contiguous, so we
    # make it so here
    if isinstance(data, np.ndarray):
        data = np.ascontiguousarray(data)
    elif isinstance(data, str):
        data = data.encode("utf-8")
    # Try to hash. Some tests return hashable types (like ndarrays) and
    # others don't (such as dictionaries)
    try:
        hd = hashlib.md5(data).hexdigest()
        # Handle those tests that return non-hashable types. This is done
        # here instead of in the tests themselves to try and reduce boilerplate
        # and provide a central location where all of this is done in case it needs
        # to be changed
    except TypeError:
        if isinstance(data, dict):
            hd = _hash_dict(data)
        elif data is None:
            hd = hashlib.md5(bytes(str(-1).encode("utf-8"))).hexdigest()
        else:
            raise
    return hd


def _save_result(data, output_file):
    r"""
    Saves the test results to the desired answer file.

    Parameters
    ----------
    data : dict
        Test results to be saved.

    output_file : str
        Name of file to save results to.
    """
    with open(output_file, "a") as f:
        yaml.dump(data, f)


def _save_raw_arrays(arrays, answer_file, func_name):
    r"""
    Saves the raw arrays produced from answer tests to a file.

    The structure of `answer_file` is: each test function (e.g.,
    test_toro1d[0-None-None-0]) forms a group. Within each group is a
    hdf5 dataset named after the test (e.g., field_values). The value
    stored in each dataset is the raw array corresponding to that
    test and function.

    Parameters
    ----------
    arrays : dict
        Keys are the test name (e.g. field_values) and the values are
        the actual answer arrays produced by the test.

    answer_file : str
        The name of the file to save the answers to, in hdf5 format.

    func_name : str
        The name of the function (possibly augmented by pytest with
        test parameters) that called the test functions
        (e.g, test_toro1d).
    """
    with h5py.File(answer_file, "a") as f:
        grp = f.create_group(func_name)
        for test_name, test_data in arrays.items():
            # Some answer tests (e.g., grid_values, projection_values)
            # return a dictionary, which cannot be handled by h5py
            if isinstance(test_data, dict):
                sub_grp = grp.create_group(test_name)
                _parse_raw_answer_dict(test_data, sub_grp)
            else:
                # Some tests return None, which hdf5 can't handle, and there is
                # no proxy, so we have to make one ourselves. Using -1
                if test_data is None:
                    test_data = -1
                grp.create_dataset(test_name, data=test_data)


def _parse_raw_answer_dict(d, h5grp):
    for k, v in d.items():
        if isinstance(v, dict):
            h5_sub_grp = h5grp.create_group(k)
            _parse_raw_answer_dict(v, h5_sub_grp)
        else:
            k = str(k)
            h5grp.create_dataset(k, data=v)


def _compare_raw_arrays(arrays, answer_file, func_name):
    r"""
    Reads in previously saved raw array data and compares the current
    results with the old ones.

    The structure of `answer_file` is: each test function (e.g.,
    test_toro1d[0-None-None-0]) forms a group. Within each group is a
    hdf5 dataset named after the test (e.g., field_values). The value
    stored in each dataset is the raw array corresponding to that
    test and function.

    Parameters
    ----------
    arrays : dict
        Keys are the test name (e.g. field_values) and the values are
        the actual answer arrays produced by the test.

    answer_file : str
        The name of the file to load the answers from, in hdf5 format.

    func_name : str
        The name of the function (possibly augmented by pytest with
        test parameters) that called the test functions
        (e.g, test_toro1d).
    """
    with h5py.File(answer_file, "r") as f:
        for test_name, new_answer in arrays.items():
            np.testing.assert_array_equal(f[func_name][test_name][:], new_answer)


def can_run_ds(ds_fn, file_check=False):
    r"""
    Validates whether or not a given input can be loaded and used as a
    Dataset object.
    """
    if isinstance(ds_fn, Dataset):
        return True
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    if file_check:
        return os.path.isfile(os.path.join(path, ds_fn))
    try:
        load(ds_fn)
        return True
    except FileNotFoundError:
        return False


def can_run_sim(sim_fn, sim_type, file_check=False):
    r"""
    Validates whether or not a given input can be used as a simulation
    time-series object.
    """
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    if file_check:
        return os.path.isfile(os.path.join(path, sim_fn))
    try:
        load_simulation(sim_fn, sim_type)
    except FileNotFoundError:
        return False
    return True


def data_dir_load(ds_fn, cls=None, args=None, kwargs=None):
    r"""
    Loads a sample dataset from the designated test_data_dir for use in
    testing.
    """
    args = args or ()
    kwargs = kwargs or {}
    path = ytcfg.get("yt", "test_data_dir")
    # Some frontends require their field_lists during test parameterization.
    # If the data isn't found, the parameterizing functions return None, since
    # pytest.skip cannot be called outside of a test or fixture.
    if ds_fn is None:
        raise FileNotFoundError
    if not os.path.isdir(path):
        raise FileNotFoundError
    if isinstance(ds_fn, Dataset):
        return ds_fn
    if cls is None:
        ds = load(ds_fn, *args, **kwargs)
    else:
        ds = cls(os.path.join(path, ds_fn), *args, **kwargs)
    ds.index
    return ds


def requires_ds(ds_fn, file_check=False):
    r"""
    Meta-wrapper for specifying required data for a test and
    checking if said data exists.
    """

    def ffalse(func):
        @functools.wraps(func)
        def skip(*args, **kwargs):
            msg = f"{ds_fn} not found, skipping {func.__name__}."
            pytest.skip(msg)

        return skip

    def ftrue(func):
        @functools.wraps(func)
        def true_wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return true_wrapper

    if not can_run_ds(ds_fn, file_check):
        return ffalse
    else:
        return ftrue


def requires_sim(sim_fn, sim_type, file_check=False):
    r"""
    Meta-wrapper for specifying a required simulation for a test and
    checking if said simulation exists.
    """

    def ffalse(func):
        @functools.wraps(func)
        def skip(*args, **kwargs):
            msg = f"{sim_fn} not found, skipping {func.__name__}."
            pytest.skip(msg)

        return skip

    def ftrue(func):
        @functools.wraps(func)
        def true_wrapper(*args, **kwargs):
            return func

        return true_wrapper

    if not can_run_sim(sim_fn, sim_type, file_check):
        return ffalse
    else:
        return ftrue


def _create_plot_window_attribute_plot(ds, plot_type, field, axis, pkwargs=None):
    r"""
    Convenience function used in plot_window_attribute_test.

    Parameters
    ----------
    ds : Dataset
        The Dataset object from which the plotting data is taken.

    plot_type : string
        Type of plot to make (e.g., SlicePlot).

    field : yt field
        The field (e.g, density) to plot.

    axis : int
        The plot axis to plot or project along.

    pkwargs : dict
        Any keywords to be passed when creating the plot.
    """
    if plot_type is None:
        raise RuntimeError("Must explicitly request a plot type")
    cls = getattr(pw, plot_type, None)
    if cls is None:
        cls = getattr(particle_plots, plot_type)
    plot = cls(ds, axis, field, **pkwargs)
    return plot


def _create_phase_plot_attribute_plot(
    data_source, x_field, y_field, z_field, plot_type, plot_kwargs=None
):
    r"""
    Convenience function used in phase_plot_attribute_test.

    Parameters
    ----------
    data_source : Dataset object
        The Dataset object from which the plotting data is taken.

    x_field : yt field
        Field to plot on x-axis.

    y_field : yt field
        Field to plot on y-axis.

    z_field : yt field
        Field to plot on z-axis.

    plot_type : string
        Type of plot to make (e.g., SlicePlot).

    plot_kwargs : dict
        Any keywords to be passed when creating the plot.
    """
    if plot_type is None:
        raise RuntimeError("Must explicitly request a plot type")
    cls = getattr(profile_plotter, plot_type, None)
    if cls is None:
        cls = getattr(particle_plots, plot_type)
    plot = cls(*(data_source, x_field, y_field, z_field), **plot_kwargs)
    return plot


def get_parameterization(fname):
    """
    Returns a dataset's field list to make test parameterizationn easier.

    Some tests (such as those that use the toro1d dataset in enzo) check
    every field in a dataset. In order to parametrize the tests without
    having to hardcode a list of every field, this function is used.
    Additionally, if the dataset cannot be found, this function enables
    pytest to mark the test as failed without the whole test run crashing,
    since the parameterization happens at import time.
    """
    try:
        ds = data_dir_load(fname)
        return ds.field_list
    except FileNotFoundError:
        return [
            None,
        ]

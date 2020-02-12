"""
Title:   utils.py
Purpose: Contains utility functions for yt answer tests
Notes:
"""
import functools
import hashlib
import inspect
import os

import numpy as np
import pytest
import yaml

from yt.config import ytcfg
from yt.convenience import load, simulation
from yt.data_objects.selection_data_containers import YTRegion
from yt.data_objects.static_output import Dataset
from yt.frontends.ytdata.api import save_as_dataset
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
from yt.utilities.exceptions import \
    YTOutputNotIdentified
import yt.visualization.particle_plots as particle_plots
import yt.visualization.plot_window as pw
import yt.visualization.profile_plotter as profile_plotter
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

    Parameters:
    -----------
        params : dict
            The dictionary of test parameters in the form
            {param_name : param_value}.

    Returns:
    --------
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
        if not isinstance(key, str) and hasattr(key, '__iter__'):
            key = _iterable_to_string(key)
        # The value can also be nested iterables
        if not isinstance(value, str) and hasattr(value, '__iter__'):
            value = _iterable_to_string(value)
        # Scene objects need special treatment to make them more IO friendly
        if isinstance(value, Scene):
            value = 'Scene' 
        elif isinstance(value, YTRegion):
            value = 'Region'
        streamlined_params[key] = value
    return streamlined_params

def _iterable_to_string(iterable):
    r"""
    An extension of streamline_for_io that does the work of making an
    iterable more io-friendly.

    Parameters:
    -----------
        iterable : python iterable
            The object to be parsed and converted.

    Returns:
    --------
        result : str
            The io-friendly version of the passed iterable.
    """
    result = iterable.__class__.__name__ 
    for elem in iterable:
        # Check for user-defined functions
        if inspect.isfunction(elem):
            result += '_' + elem.__name__
        # Non-string iterables (e.g., lists, tuples, etc.)
        elif not isinstance(elem, str) and hasattr(elem, '__iter__'):
            result += '_' + _iterable_to_string(elem)
        # Non-string non-iterables (ints, floats, etc.)
        elif not isinstance(elem, str) and not hasattr(elem, '__iter__'):
            result += '_' + str(elem)
        # Strings
        elif isinstance(elem, str):
            result += '_' + elem
    return result

def _hash_results(results):
    r"""
    Driver function for hashing the test result.

    Parameters:
    -----------
        results : dict
            Dictionary of {test_name : test_result} pairs.

    Returns:
    --------
        results : dict
            Same as the passed results, but the test_results are now
            hex digests of the hashed test_result.
    """
    # Here, results should be comprised of only the tests, not the test
    # parameters
    for test_name, test_value in results.items():
        # These tests have issues with python-specific anchors and so
        # are already hashed
        # (see their definitions in yt/utilites/answer_testing/answer_tests.py)
        if test_name in ['projection_values', 'pixelized_projection_values', 'grid_values']:
            continue
        else:
            results[test_name] = generate_hash(test_value)
    return results

def _hash_dict(data):
    r"""
    Specifically handles hashing a dictionary object. 

    Parameters:
    -----------
        data : dict
            The dictionary to be hashed.

    Returns:
    --------
        hd.hexdigest : str
            The hex digest of the hashed dictionary.
    """
    hd = None
    for key, value in sorted(data.items()):
        if hd is None:
            hd = hashlib.md5(bytearray(key.encode('utf8')) + bytearray(value))
        else:
            hd.update(bytearray(key.encode('utf8')) + bytearray(value))
    return hd.hexdigest()

def generate_hash(data):
    r"""
    Actually performs the hash operation. 

    Parameters:
    -----------
        data : python object 
            Data to be hashed.

    Returns:
    --------
        hd : str
            Hex digest of the hashed data.
    """
    # Sometimes md5 complains that the data is not contiguous, so we
    # make it so here
    if isinstance(data, np.ndarray):
        data = np.ascontiguousarray(data)
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
        else:
            raise TypeError
    return hd

def _save_result(data, outputFile):
    r"""
    Saves the test results to the desired answer file. 

    Parameters:
    -----------
        data : dict
            Test results to be saved.

        outputFile : str
            Name of file to save results to.

    Returns:
    --------
        None 
    """
    with open(outputFile, 'a') as f:
        yaml.dump(data, f)

def _compare_result(data, outputFile):
    r"""
    Compares the just-generated test results to those that are already
    saved.

    Parameters:
    -----------
        data : dict
            Just-generated test results.

        outputFile : str
            Name of file where answers are already saved.
    """
    # Load the saved data
    with open(outputFile, 'r') as f:
        savedData = yaml.safe_load(f)
    # Define the comparison function
    def _check_vals(newVals, oldVals):
        for key, value in newVals.items():
            if isinstance(value, dict):
                _check_vals(value, oldVals[key])
            else:
                assert value == oldVals[key]
    # Compare
    _check_vals(data, savedData)

def _handle_hashes(save_dir_name, fname, hashes, answer_store):
    r"""
    Driver function for deciding whether to save the test results or
    compare them to already saved results.

    Parameters:
    -----------
        save_dir_name : str
            Name of the directory to save results or where results are
            already saved.

        fname : str
            Name of the file to either save results to or where results
            are already saved.

        hashes : dict
            The just-generated test results.

        answer_store : bool
            If true, save the just-generated test results, otherwise,
            compare them to the previously saved results.
    """
    # Set up the answer file in the answer directory
    answer_file = os.path.join(save_dir_name, fname)
    # Save the result
    if answer_store:
        _save_result(hashes, answer_file)
    # Compare to already saved results
    else:
        _compare_result(hashes, answer_file)


def _save_arrays(save_dir_name, fbasename, arrays, answer_store):
    r"""
    Driver routine for either saving the raw arrays resulting from the
    tests, or compare them to previously saved results.

    Parameters:
    -----------
        save_dir_name : str
            Name of the directory to save results or where results are
            already saved.

        fbasename : str
            Base name (no extension) of the file to either save results
            to or where results are already saved.

        arrays : dict
            The raw arrays generated from the tests, with the test name
            as a key.

        answer_store : bool
            If true, save the just-generated test results, otherwise,
            compare them to the previously saved results.
    """
    pass


def can_run_ds(ds_fn, file_check = False):
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
    except YTOutputNotIdentified:
        return False


def can_run_sim(sim_fn, sim_type, file_check = False):
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
        simulation(sim_fn, sim_type)
    except YTOutputNotIdentified:
        return False
    return True 


def data_dir_load(ds_fn, cls = None, args = None, kwargs = None):
    r"""
    Loads a sample dataset from the designated test_data_dir for use in
    testing.
    """
    args = args or ()
    kwargs = kwargs or {}
    path = ytcfg.get("yt", "test_data_dir")
    if isinstance(ds_fn, Dataset): return ds_fn
    if not os.path.isdir(path):
        return False
    if cls is None:
        ds = load(ds_fn, *args, **kwargs)
    else:
        ds = cls(os.path.join(path, ds_fn), *args, **kwargs)
    ds.index
    return ds


def requires_ds(ds_fn, file_check = False):
    r"""
    Meta-wrapper for specifying required data for a test and
    checking if said data exists.
    """
    def ffalse(func):
        @functools.wraps(func)
        def skip(*args, **kwargs):
            msg = "{} not found, skipping {}.".format(ds_fn, func.__name__)
            pytest.fail(msg)
        return skip
    def ftrue(func):
        @functools.wraps(func)
        def true_wrapper(*args, **kwargs):
            return func
        return true_wrapper
    if not can_run_ds(ds_fn, file_check):
        return ffalse
    else:
        return ftrue


def requires_sim(sim_fn, sim_type, file_check = False):
    r"""
    Meta-wrapper for specifying a required simulation for a test and
    checking if said simulation exists.
    """
    def ffalse(func):
        @functools.wraps(func)
        def skip(*args, **kwargs):
            msg = "{} not found, skipping {}.".format(sim_fn, func.__name__)
            pytest.fail(msg)
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


def create_obj(ds, obj_type):
    # obj_type should be tuple of
    #  ( obj_name, ( args ) )
    if obj_type is None:
        return ds.all_data()
    cls = getattr(ds, obj_type[0])
    obj = cls(*obj_type[1])
    return obj


def compare_unit_attributes(ds1, ds2):
    r"""
    Checks to make sure that the length, mass, time, velocity, and
    magnetic units are the same for two different dataset objects.
    """
    attrs = ('length_unit', 'mass_unit', 'time_unit',
             'velocity_unit', 'magnetic_unit')
    for attr in attrs:
        u1 = getattr(ds1, attr, None)
        u2 = getattr(ds2, attr, None)
        assert u1 == u2


def fake_halo_catalog(data):
    filename = "catalog.0.h5"
    ftypes = dict((field, '.') for field in data)
    extra_attrs = {"data_type": "halo_catalog",
                   "num_halos": data['particle_mass'].size}
    ds = {'cosmological_simulation': 1,
          'omega_lambda': 0.7,
          'omega_matter': 0.3,
          'hubble_constant': 0.7,
          'current_redshift': 0,
          'current_time': YTQuantity(1, 'yr'),
          'domain_left_edge': YTArray(np.zeros(3), 'cm'),
          'domain_right_edge': YTArray(np.ones(3), 'cm')}
    save_as_dataset(ds, filename, data, field_types=ftypes,
        extra_attrs=extra_attrs
    )
    return filename


def _create_plot_window_attribute_plot(ds, ptype, field, axis, pkwargs = None):
    r"""
    Convenience function used in plot_window_attribute_test.

    Parameters:
    -----------
        ds : Dataset
            The Dataset object from which the plotting data is taken.

        ptype : string
            Type of plot to make (e.g., SlicePlot).

        field : yt field
            The field (e.g, density) to plot.

        axis : int
            The plot axis to plot or project along.
            
        pkwargs : dict
            Any keywords to be passed when creating the plot.
    """
    if ptype is None:
        raise RuntimeError('Must explicitly request a plot type')
    cls = getattr(pw, ptype, None)
    if cls is None:
        cls = getattr(particle_plots, ptype)
    plot = cls(*(ds, axis, field), **pkwargs)
    return plot


def _create_phase_plot_attribute_plot(data_source, x_field, y_field, z_field,
                plot_type, plot_kwargs=None):
    r"""
    Convenience function used in phase_plot_attribute_test.

    Parameters:
    -----------
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
        raise RuntimeError('Must explicitly request a plot type')
    cls = getattr(profile_plotter, plot_type, None)
    if cls is None:
        cls = getattr(particle_plots, plot_type)
    plot = cls(*(data_source, x_field, y_field, z_field), **plot_kwargs)
    return plot

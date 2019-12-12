"""
Title:   utils.py
Purpose: Contains utility functions for yt answer tests
Notes:
"""
import hashlib
import os

import numpy as np
import pytest
import yaml

from yt.config import ytcfg
from yt.convenience import load, simulation
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


#============================================
#            streamline_for_io
#============================================
def streamline_for_io(params):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    streamlined_params = {}
    for key, value in params.items():
        # The key can be nested iterables, e.g., 
        # d = [None, ('sphere', (center, (0.1, 'unitary')))] so we need
        # to use recursion
        if not isinstance(key, str) and hasattr(key, '__iter__'):
            key = iterable_to_string(key)
        # The value can also be nested iterables
        if not isinstance(value, str) and hasattr(value, '__iter__'):
            value = iterable_to_string(value)
        streamlined_params[key] = value
    return streamlined_params

#============================================
#            iterable_to_string
#============================================
def iterable_to_string(iterable):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    result = iterable.__class__.__name__ 
    for elem in iterable:
        # Non-string iterables (e.g., lists, tuples, etc.)
        if not isinstance(elem, str) and hasattr(elem, '__iter__'):
            result += '_' + iterable_to_string(elem)
        # Non-string non-iterables (ints, floats, etc.)
        elif not isinstance(elem, str) and not hasattr(elem, '__iter__'):
            result += '_' + str(elem)
        # Strings
        elif isinstance(elem, str):
            result += '_' + elem
    return result

#============================================
#               hash_results
#============================================
def hash_results(results):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    # Here, results should be comprised of only the tests, not the test
    # parameters
    for test_name, test_value in results.items():
        # These tests have issues with python-specific anchors and so
        # are already hashed
        # (see their definitions in yt/utilites/answer_testing/framework.py)
        if test_name in ['projection_values', 'pixelized_projection_values', 'grid_values']:
            continue
        else:
            results[test_name] = generate_hash(test_value)
    return results

#============================================
#                hash_dict
#============================================
def hash_dict(data):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    hd = None
    for key, value in data.items():
        if hd is None:
            hd = hashlib.md5(bytes(key.encode('utf8')) + value.tobytes())
        else:
            hd.update(bytes(key.encode('utf8')) + value.tobytes())
    return hd.hexdigest()

#============================================
#              generate_hash
#============================================
def generate_hash(data):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    if isinstance(data, np.ndarray):
        # Sometimes md5 complains that the data is not contiguous
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
            hd = hash_dict(data)
        else:
            raise TypeError
    return hd

#============================================
#                 save_result
#============================================
def save_result(data, outputFile):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    with open(outputFile, 'a') as f:
        yaml.dump(data, f)

#============================================
#               compare_result
#============================================
def compare_result(data, outputFile):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    # Load the saved data
    with open(outputFile, 'r') as f:
        savedData = yaml.safe_load(f)
    # Define the comparison function
    def check_vals(newVals, oldVals):
        for key, value in newVals.items():
            if isinstance(value, dict):
                check_vals(value, oldVals[key])
            else:
                assert value == oldVals[key]
    # Compare
    check_vals(data, savedData)

#============================================
#               handle_hashes
#============================================
def handle_hashes(save_dir_name, fname, hashes, answer_store):
    """
    Doc string.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    # Set up the answer file in the answer directory
    answer_file = os.path.join(save_dir_name, fname)
    # Save the result
    if answer_store:
        save_result(hashes, answer_file)
    # Compare to already saved results
    else:
        compare_result(hashes, answer_file)


#============================================
#                 can_run_ds
#============================================
def can_run_ds(ds_fn, file_check = False):
    """
    Checks to see if the given file name is a valid data set.

    Parameters:
    -----------
        pass

    Raises:
    -------
        returns
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

#============================================
#                can_run_sim
#============================================
def can_run_sim(sim_fn, sim_type, file_check = False):
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

#============================================
#                data_dir_load
#============================================
def data_dir_load(ds_fn, cls = None, args = None, kwargs = None):
    """
    Loads the given data set.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
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

#============================================-
#                 requires_ds
#============================================
def requires_ds(ds_fn, file_check = False):
    """
    Meta-wrapper for specifying required data for a test and
    checking if said data exists.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    def ffalse(func):
        def skip(*args, **kwargs):
            msg = "{} not found, skipping {}.".format(ds_fn, func.__name__)
            pytest.fail(msg)
        return skip
    def ftrue(func):
        return func
    if not can_run_ds(ds_fn, file_check):
        return ffalse
    else:
        return ftrue

#============================================
# requires_sim
#============================================
def requires_sim(sim_fn, sim_type, file_check = False):
    def ffalse(func):
        def skip(*args, **kwargs):
            msg = "{} not found, skipping {}.".format(sim_fn, func.__name__)
            pytest.fail(msg)
        return skip
    def ftrue(func):
        return func
    if not can_run_sim(sim_fn, sim_type, file_check):
        return ffalse
    else:
        return ftrue


#============================================
#                create_obj
#============================================
def create_obj(ds, obj_type):
    """
    Builds a dataset object of the desired type.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    # obj_type should be tuple of
    #  ( obj_name, ( args ) )
    if obj_type is None:
        return ds.all_data()
    cls = getattr(ds, obj_type[0])
    obj = cls(*obj_type[1])
    return obj


#============================================
#          compare_unit_attributes
#============================================
def compare_unit_attributes(ds1, ds2):
    """
    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    attrs = ('length_unit', 'mass_unit', 'time_unit',
             'velocity_unit', 'magnetic_unit')
    for attr in attrs:
        u1 = getattr(ds1, attr, None)
        u2 = getattr(ds2, attr, None)
        assert u1 == u2


#============================================
#              fake_halo_catalog
#============================================
def fake_halo_catalog(data):
    """
    Helper function to create and save a mock halo catalog for use
    in testing.

    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
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


#============================================
#                 create_plot
#+===========================================
def create_plot(ds, plot_type, plot_field, plot_axis, plot_kwargs = None):
    """
    Parameters:
    -----------
        pass

    Raises:
    -------
        pass

    Returns:
    --------
        pass
    """
    # plot_type should be a string
    # plot_kwargs should be a dict
    if plot_type is None:
        raise RuntimeError('Must explicitly request a plot type')
    cls = getattr(pw, plot_type, None)
    if cls is None:
        cls = getattr(particle_plots, plot_type)
    plot = cls(*(ds, plot_axis, plot_field), **plot_kwargs)
    return plot


#============================================
#               create_plot2
#============================================
def create_plot2(data_source, x_field, y_field, z_field,
                plot_type, plot_kwargs=None):
    # plot_type should be a string
    # plot_kwargs should be a dict
    if plot_type is None:
        raise RuntimeError('Must explicitly request a plot type')
    cls = getattr(profile_plotter, plot_type, None)
    if cls is None:
        cls = getattr(particle_plots, plot_type)
    plot = cls(*(data_source, x_field, y_field, z_field), **plot_kwargs)
    return plot

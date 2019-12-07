"""
Title:   utils.py
Purpose: Contains utility functions for yt answer tests
Notes:
"""
from collections import OrderedDict
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
#               array_to_hash
#============================================
def array_to_hash(d):
    """
    This function loops recursively over each nested dictionary in d
    and, when it reaches a non-dictionary value, if that value is an
    array, it converts it to a bytes array.
    """
    for k, v in d.items():
        if k == 'grid_values' or k == 'pixelized_projection_values':
            continue
        if isinstance(v, dict) or isinstance(v, OrderedDict):
            array_to_hash(v)
        else:
            d[k] = generate_hash(v)
            if not isinstance(k, str):
                d.update({k.__repr__() : d[k]})
                del d[k]
    return d


#============================================
#               generate_hash
#============================================
def generate_hash(data):
    """
    Generates an md5 hash from the given data. This function assumes
    that the data passed with either be a dataset object or else the
    result of a particular test. It assumes that the test result has
    already been put into a hashable form, since the data returned by
    each test is different.

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
    try:
        hd = hashlib.md5(data).hexdigest()
    except TypeError:
        hd = data.__repr__()
    return hd


#============================================
#               log_test_error
#============================================
def log_test_error(func_name, val, param, saved_hashes, e):
    """
    Logs the answer test error by saving the relevant information for
    both the saved and obtained tests.

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
    msg = "Error: {} in {} with value {} differs from saved value of {}".format(
        param, func_name, val, saved_hashes[func_name][param])
    pytest.fail(msg)


#============================================
#                 check_vals
#============================================
def check_vals(newVals, savedVals):
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
    for key, value in newVals.items():
        if isinstance(value, dict):
            check_vals(value, savedVals[key])
        else:
            assert value == savedVals[key]


#============================================
#               handle_hashes
#============================================
def handle_hashes(save_dir_name, fname, hashes, answer_store):
    """
    Either saves the answers for later comparison or loads in the saved
    answers and does the comparison.

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
    fname = os.path.join(save_dir_name, fname)
    # Save answer
    if answer_store:
        with open(fname, 'a') as f:
            yaml.dump(hashes, f, default_flow_style=False, Dumper=yaml.Dumper)
    # Compare to already saved answer
    else:
        with open(fname, 'r') as f:
            saved_hashes = yaml.load(f, Loader=yaml.Loader)
        # The layout of the answer file is:
        # {functionName : {param1 : val1, param2 : val2, ...}, functionName2 : {...}}
        # Since dicts are being used we can't bank on the test parameters being in
        # the same order in the file and the hashes dict. OrderedDicts are not used
        # because there's a bug in pyyaml with duplicate anchors. It's also harder
        # to read
        import pdb; pdb.set_trace()
        check_vals(hashes, saved_hashes)
        # for func_name, test_vals in hashes.items():
        #     for test_param, val in test_vals.items():
        #         try:
        #             assert val == saved_hashes[func_name][test_param]
        #         except (AssertionError, KeyError) as e:
        #             log_test_error(func_name, val, test_param, saved_hashes, e)


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
            print("{} not found, skipping {}.".format(ds_fn, func.__name__))
            assert False
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
            print("{} not found, skipping {}.".format(sim_fn, func.__name__))
            assert False
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

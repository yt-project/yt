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
from yt.convenience import load
from yt.data_objects.static_output import Dataset
from yt.frontends.ytdata.api import save_as_dataset
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
from yt.utilities.exceptions import \
    YTOutputNotIdentified
import yt.visualization.particle_plots as particle_plots
import yt.visualization.profile_plotter as profile_plotter


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
    return hashlib.md5(data).hexdigest()


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
            yaml.dump(hashes, f, default_flow_style=False)
    # Compare to already saved answer
    else:
        with open(fname, 'r') as f:
            saved_hashes = yaml.load(fname)
        assert hashes == saved_hashes


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
        return lambda: None
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
        return lambda: None
    def ftrue(func):
        return func
    elif not can_run_sim(sim_fn, sim_type, file_check):
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
#         convert_to_ordered_dict
#============================================
def convert_to_ordered_dict(d):
    """
    Converts the given dictionary to an OrderedDict sorted 
    alphabetically by key.

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
    # First store the key-value pairs in a list and sort by key. When
    # the entries of the list being sorted are tuples, sorted()
    # defaults to using the first entry of the tuple. This trend goes
    # on if the tuples are nested (e.g., if the dict is:
    # d = {('enzo', 'velocity') : val2, ('enzo', 'density') : val1, }
    # then the below code will give: l = [(('enzo', 'density'), val1),
    # (('enzo', 'velocity'), val2)]
    l = sorted([(k, v) for k, v in d.items()])
    # Put into an OrderedDict
    od = OrderedDict()
    for entry in l:
        od[entry[0]] = entry[1]
    return od


#============================================
#          check_result_hashability
#============================================
def check_result_hashability(result):
    """
    This is specifically for the parentage_relationships test.  
    There are three cases: element 0 is empty and element 1 isn't, 
    vice versa, or neither is empty, but they don't have the same
    length. For whatever reason, the hexes for a np.array hash are
    only the same if the rows in the array have the same length.
    That is, a = np.array([[1,2,3], [4]]) and
    b = np.array([[1,2,3], [4]]) will have different hashes
    despite containing the same data. This holds for all
    configurations unless the rows of the array have the same
    number of elements (even both empty is fine). This pads the array
    that's too short with -2, since I used -1 for None in the actual
    test. This function generalizes from the two grids used in kh2d to
    an arbitrary number of grids. This only applies to the children
    key, because result["children"] = [list1, list2]. For
    result["parents"], that's just a list of ints corresponding to grid
    ids, so it doesn't need to be changed.

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
    # Get longest element
    longest = 0
    for sublist in result["children"]:
        if len(sublist) > longest:
            longest = len(sublist)
    # Now adjust such that each sublist has the same length as the
    # longest one
    for sublist in result["children"]:
        if len(sublist) < longest:
            diff = longest - len(sublist)
            appendix = [-2 for i in range(diff)]
            sublist += appendix
    return result


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
def create_plot(self, ds, plot_type, plot_field, plot_axis, plot_kwargs = None):
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
def create_plot(data_source, x_field, y_field, z_field,
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

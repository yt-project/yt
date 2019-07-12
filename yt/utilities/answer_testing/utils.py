"""
Title:   utils.py
Purpose: Contains utility functions for yt answer tests
Notes:
"""
from collections import OrderedDict
import hashlib
import os

from yt.config import ytcfg
from yt.convenience import load
from yt.data_objects.static_output import Dataset
from yt.utilities.exceptions import \
    YTOutputNotIdentified


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
#               store_hashes
#============================================
def store_hashes(ds_name, hashes):
    """
    This writes the hashes for the data contained in the data file
    and the hashes for the results of each test run on that data
    file.

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
    with open(ds_name + '-hashes.txt', 'w') as f:
        # Write the hashes to the file
        for key, value in hashes.items():
            f.write(key + ' : ' + str(value) + '\n')

#============================================
#                load_hashes
#============================================
def load_hashes(ds_name):
    """
    Reads the hashes for each test related to a given data set and
    returns them.

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
    hashes = {}
    with open(ds_name + '-hashes.txt', 'r') as f:
        for line in f:
            key, value = line.split(':')
            hashes[key.strip()] = value.strip()
    return hashes


#============================================
#               handle_hashes
#============================================
def handle_hashes(ds_name, hashes, answer_store):
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
    # Save answer
    if answer_store:
        store_hashes(ds_name, hashes)
    # Compare to already saved answer
    else:
        saved_hashes = load_hashes(ds_name)
        for key in hashes.keys():
            try:
                assert saved_hashes[key] == hashes[key]
            except AssertionError:
                print("Error, hashes for ds: %s test: %s don't match!" %
                    (ds_name, key)
                )


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
def requires_ds(ds_fn, big_data = False, file_check = False):
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
    run_big_data = False
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if run_big_data is False and big_data is True:
        return ffalse
    elif not can_run_ds(ds_fn, file_check):
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

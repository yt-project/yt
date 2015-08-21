"""
Utility functions for ytdata frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np

from yt.funcs import \
    mylog
from yt.units.yt_array import \
    YTArray

def to_yt_dataset(ds, filename, data, field_types=None,
                  extra_attrs=None):
    r"""Export a set of field arrays to a reloadable yt dataset.

    This function can be used to create a yt loadable dataset from a 
    set of arrays.  The field arrays can either be associated with a 
    loaded dataset or, if not, a dictionary of dataset attributes can
    be provided that will be used as metadata for the new dataset.  The 
    resulting dataset can be reloaded as a yt dataset.

    Parameters
    ----------
    ds : dataset
        The dataset associated with the fields.  
    filename : str
        The name of the file to be written.
    data : dict
        A dictionary of field arrays to be saved.
    extra_attrs: dict
        A dictionary of additional attributes to be saved.

    Returns
    -------
    filename : str
        The name of the file that has been created.

    Examples
    --------

    >>> import yt
    >>> from yt.frontends.ytdata.api import to_yt_dataset
    >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> sphere = ds.sphere([0.5]*3, (10, "Mpc")
    >>> sphere_density = sphere["density"]
    >>> region = ds.box([0.]*3, [0.25]*3)
    >>> region_density = region["density"]
    >>> data = {}
    >>> data["sphere_density"] = sphere_density
    >>> data["region_density"] = region_density
    >>> to_yt_dataset(ds, "density_data.h5", data)

    >>> import yt
    >>> from yt.frontends.ytdata.api import to_yt_dataset
    >>> from yt.units.yt_array import YTArray, YTQuantity
    >>> data = {"density": YTArray(np.random.random(10), "g/cm**3"),
    ...         "temperature": YTArray(np.random.random(10), "K")}
    >>> ds_data = {"domain_left_edge": YTArray(np.zeros(3), "cm"),
    ...            "domain_right_edge": YTArray(np.ones(3), "cm"),
    ...            "current_time": YTQuantity(10, "Myr")}
    >>> to_yt_dataset(ds_data, "random_data.h5", data)
    
    """

    mylog.info("Saving field data to yt dataset: %s." % filename)

    if extra_attrs is None: extra_attrs = {}
    base_attrs  = ["domain_left_edge", "domain_right_edge",
                   "current_redshift", "current_time",
                   "domain_dimensions", "periodicity",
                   "cosmological_simulation", "omega_lambda",
                   "omega_matter", "hubble_constant"]

    fh = h5py.File(filename, "w")
    for attr in base_attrs:
        if isinstance(ds, dict):
            my_val = ds.get(attr, None)
        else:
            my_val = getattr(ds, attr, None)
        if my_val is None:
            mylog.warn("Skipping %s attribute, this may be just fine." % attr)
            continue
        if hasattr(my_val, "units"):
            my_val = my_val.in_cgs()
        fh.attrs[attr] = my_val
    for attr in extra_attrs:
        my_val = extra_attrs[attr]
        if hasattr(my_val, "units"):
            my_val = my_val.in_cgs()
        fh.attrs[attr] = my_val
    if "data_type" not in extra_attrs:
        fh.attrs["data_type"] = "yt_array_data"
    for field in data:
        if field_types is None:
            field_type = "data"
        else:
            field_type = field_types[field]
        if field_type not in fh:
            fh.create_group(field_type)
        
        # for now, let's avoid writing "code" units
        if hasattr(field, "units"):
            data[field].convert_to_cgs()
        dataset = _yt_array_hdf5(fh[field_type], field, data[field])
    fh.close()

def _hdf5_yt_array(fh, field, ds=None):
    r"""Load an hdf5 dataset as a YTArray.

    Reads in a dataset from an open hdf5 file or group and uses the
    "units" attribute, if it exists, to apply units.
    
    Parameters
    ----------
    fh : an open hdf5 file or hdf5 group
        The hdf5 file or group in which the dataset exists.
    field : str
        The name of the field to be loaded.
    ds : yt Dataset
        If not None, the unit_registry of the dataset
        is used to apply units.

    Returns
    -------
    A YTArray of the requested field.
    
    """
    
    if ds is None:
        new_arr = YTArray
    else:
        new_arr = ds.arr
    units = ""
    if "units" in fh[field].attrs:
        units = fh[field].attrs["units"]
    if units == "dimensionless": units = ""
    return new_arr(fh[field].value, units)

def _yt_array_hdf5(fh, field, data):
    r"""Save a YTArray to an open hdf5 file or group.

    Save a YTArray to an open hdf5 file or group, and save the 
    units to a "units" attribute.
    
    Parameters
    ----------
    fh : an open hdf5 file or hdf5 group
        The hdf5 file or group to which the data will be written.
    field : str
        The name of the field to be saved.
    ddata : YTArray
        The data array to be saved.

    Returns
    -------
    dataset : hdf5 dataset
        The created hdf5 dataset.
    
    """

    dataset = fh.create_dataset(str(field), data=data)
    units = ""
    if isinstance(data, YTArray):
        units = str(data.units)
    dataset.attrs["units"] = units
    return dataset

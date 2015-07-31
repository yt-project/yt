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

from yt.funcs import mylog

def to_yt_dataset(ds, filename, data):
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

    Returns
    -------
    filename : str
        The name of the file that has been created.

    Examples
    --------

    COMING SOON!
    
    """

    fh = h5py.file(filename, "w")

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

    dataset = fh.create_dataset(field, data=data)
    units = ""
    if isinstance(data, YTArray):
        units = str(data.units)
    dataset.attrs["units"] = units
    return dataset

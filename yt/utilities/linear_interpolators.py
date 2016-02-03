"""
A collection of helper functions, most generally for things
that SciPy doesn't have that I expected it to



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
import yt.utilities.lib.interpolators as lib

class UnilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names, truncate=False):
        r"""Initialize a 1D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple or array
            If a tuple, this should specify the upper and lower bounds 
            for the bins of the data table.  This assumes the bins are 
            evenly spaced.  If an array, this specifies the bins 
            explicitly.
        field_names: str
            Name of the field to be used as input data for interpolation.
        truncate : bool
            If False, an exception is raised if the input values are 
            outside the bounds of the table.  If True, extrapolation is 
            performed.
        
        Examples
        --------

        ad = ds.all_data()
        table_data = np.random.random(64)
        interp = UnilinearFieldInterpolator(table_data, (0.0, 1.0), "x",
                                            truncate=True)
        field_data = interp(ad)
        
        """
        self.table = table.astype('float64')
        self.truncate = truncate
        self.x_name = field_names
        if isinstance(boundaries, np.ndarray):
            if boundaries.size != table.shape[0]:
                mylog.error("Bins array not the same length as the data.")
                raise ValueError
            self.x_bins = boundaries
        else:
            x0, x1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0]).astype('float64')
        
    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel().astype('float64')

        x_i = (np.digitize(x_vals, self.x_bins) - 1).astype('int32')
        if np.any((x_i == -1) | (x_i == len(self.x_bins)-1)):
            if not self.truncate:
                mylog.error("Sorry, but your values are outside" + \
                            " the table!  Dunno what to do, so dying.")
                mylog.error("Error was in: %s", data_object)
                raise ValueError
            else:
                x_i = np.minimum(np.maximum(x_i,0), len(self.x_bins)-2)

        my_vals = np.zeros(x_vals.shape, dtype='float64')
        lib.UnilinearlyInterpolate(self.table, x_vals, self.x_bins, x_i, my_vals)
        my_vals.shape = orig_shape
        return my_vals

class BilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names, truncate=False):
        r"""Initialize a 2D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple
            Either a tuple of lower and upper bounds for the x and y bins 
            given as (x0, x1, y0, y1) or a tuple of two arrays containing the 
            x and y bins.
        field_names: list
            Names of the fields to be used as input data for interpolation.
        truncate : bool
            If False, an exception is raised if the input values are 
            outside the bounds of the table.  If True, extrapolation is 
            performed.
        
        Examples
        --------

        ad = ds.all_data()
        table_data = np.random.random((64, 64))
        interp = BilinearFieldInterpolator(table_data, (0.0, 1.0, 0.0, 1.0), 
                                           ["x", "y"],
                                           truncate=True)
        field_data = interp(ad)
        
        """
        self.table = table.astype('float64')
        self.truncate = truncate
        self.x_name, self.y_name = field_names
        if len(boundaries) == 4:
            x0, x1, y0, y1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0]).astype('float64')
            self.y_bins = np.linspace(y0, y1, table.shape[1]).astype('float64')
        elif len(boundaries) == 2:
            if boundaries[0].size != table.shape[0]:
                mylog.error("X bins array not the same length as the data.")
                raise ValueError
            if boundaries[1].size != table.shape[1]:
                mylog.error("Y bins array not the same length as the data.")
                raise ValueError
            self.x_bins = boundaries[0]
            self.y_bins = boundaries[1]
        else:
            mylog.error("Boundaries must be given as (x0, x1, y0, y1) or as (x_bins, y_bins)")
            raise ValueError

    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel().astype('float64')
        y_vals = data_object[self.y_name].ravel().astype('float64')

        x_i = (np.digitize(x_vals, self.x_bins) - 1).astype('int32')
        y_i = (np.digitize(y_vals, self.y_bins) - 1).astype('int32')
        if np.any((x_i == -1) | (x_i == len(self.x_bins)-1)) \
            or np.any((y_i == -1) | (y_i == len(self.y_bins)-1)):
            if not self.truncate:
                mylog.error("Sorry, but your values are outside" + \
                            " the table!  Dunno what to do, so dying.")
                mylog.error("Error was in: %s", data_object)
                raise ValueError
            else:
                x_i = np.minimum(np.maximum(x_i,0), len(self.x_bins)-2)
                y_i = np.minimum(np.maximum(y_i,0), len(self.y_bins)-2)

        my_vals = np.zeros(x_vals.shape, dtype='float64')
        lib.BilinearlyInterpolate(self.table,
                                 x_vals, y_vals, self.x_bins, self.y_bins,
                                 x_i, y_i, my_vals)
        my_vals.shape = orig_shape
        return my_vals

class TrilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names, truncate = False):
        r"""Initialize a 3D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple
            Either a tuple of lower and upper bounds for the x, y, and z bins 
            given as (x0, x1, y0, y1, z0, z1) or a tuple of three arrays 
            containing the x, y, and z bins.
        field_names: list
            Names of the fields to be used as input data for interpolation.
        truncate : bool
            If False, an exception is raised if the input values are 
            outside the bounds of the table.  If True, extrapolation is 
            performed.
        
        Examples
        --------

        ad = ds.all_data()
        table_data = np.random.random((64, 64, 64))
        interp = BilinearFieldInterpolator(table_data, 
                                           (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), 
                                           ["x", "y", "z"],
                                           truncate=True)
        field_data = interp(ad)
        
        """
        self.table = table.astype('float64')
        self.truncate = truncate
        self.x_name, self.y_name, self.z_name = field_names
        if len(boundaries) == 6:
            x0, x1, y0, y1, z0, z1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0]).astype('float64')
            self.y_bins = np.linspace(y0, y1, table.shape[1]).astype('float64')
            self.z_bins = np.linspace(z0, z1, table.shape[2]).astype('float64')
        elif len(boundaries) == 3:
            if boundaries[0].size != table.shape[0]:
                mylog.error("X bins array not the same length as the data.")
                raise ValueError
            if boundaries[1].size != table.shape[1]:
                mylog.error("Y bins array not the same length as the data.")
                raise ValueError
            if boundaries[2].size != table.shape[2]:
                mylog.error("Z bins array not the same length as the data.")
                raise ValueError
            self.x_bins = boundaries[0]
            self.y_bins = boundaries[1]
            self.z_bins = boundaries[2]
        else:
            mylog.error("Boundaries must be given as (x0, x1, y0, y1, z0, z1) or as (x_bins, y_bins, z_bins)")
            raise ValueError
        
    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel().astype('float64')
        y_vals = data_object[self.y_name].ravel().astype('float64')
        z_vals = data_object[self.z_name].ravel().astype('float64')

        x_i = np.digitize(x_vals, self.x_bins).astype("int") - 1
        y_i = np.digitize(y_vals, self.y_bins).astype("int") - 1
        z_i = np.digitize(z_vals, self.z_bins).astype("int") - 1
        if np.any((x_i == -1) | (x_i == len(self.x_bins)-1)) \
            or np.any((y_i == -1) | (y_i == len(self.y_bins)-1)) \
            or np.any((z_i == -1) | (z_i == len(self.z_bins)-1)):
            if not self.truncate:
                mylog.error("Sorry, but your values are outside" + \
                            " the table!  Dunno what to do, so dying.")
                mylog.error("Error was in: %s", data_object)
                raise ValueError
            else:
                x_i = np.minimum(np.maximum(x_i,0), len(self.x_bins)-2)
                y_i = np.minimum(np.maximum(y_i,0), len(self.y_bins)-2)
                z_i = np.minimum(np.maximum(z_i,0), len(self.z_bins)-2)

        my_vals = np.zeros(x_vals.shape, dtype='float64')
        lib.TrilinearlyInterpolate(self.table,
                                 x_vals, y_vals, z_vals,
                                 self.x_bins, self.y_bins, self.z_bins,
                                 x_i, y_i, z_i, my_vals)
        my_vals.shape = orig_shape
        return my_vals

def get_centers(ds, filename, center_cols, radius_col, unit='1'):
    """
    Return an iterator over EnzoSphere objects generated from the appropriate 
    columns in *filename*.  Optionally specify the *unit* radius is in.
    """
    for line in open(filename):
        if line.startswith("#"): continue
        vals = line.split()
        x,y,z = [float(vals[i]) for i in center_cols]
        r = float(vals[radius_col])
        yield ds.sphere([x,y,z], r/ds[unit])

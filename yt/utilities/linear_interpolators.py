"""
A collection of helper functions, most generally for things
that SciPy doesn't have that I expected it to

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

from yt.funcs import *
import yt.utilities.lib as lib

class UnilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names, truncate=False):
        self.table = table.astype('float64')
        self.truncate = truncate
        x0, x1 = boundaries
        self.x_name = field_names
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
        self.table = table.astype('float64')
        self.truncate = truncate
        x0, x1, y0, y1 = boundaries
        self.x_name, self.y_name = field_names
        self.x_bins = np.linspace(x0, x1, table.shape[0]).astype('float64')
        self.y_bins = np.linspace(y0, y1, table.shape[1]).astype('float64')

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
        self.table = table.astype('float64')
        self.truncate = truncate
        x0, x1, y0, y1, z0, z1 = boundaries
        self.x_name, self.y_name, self.z_name = field_names
        self.x_bins = np.linspace(x0, x1, table.shape[0]).astype('float64')
        self.y_bins = np.linspace(y0, y1, table.shape[1]).astype('float64')
        self.z_bins = np.linspace(z0, z1, table.shape[2]).astype('float64')

    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel().astype('float64')
        y_vals = data_object[self.y_name].ravel().astype('float64')
        z_vals = data_object[self.z_name].ravel().astype('float64')

        x_i = np.digitize(x_vals, self.x_bins) - 1
        y_i = np.digitize(y_vals, self.y_bins) - 1
        z_i = np.digitize(z_vals, self.z_bins) - 1
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

def get_centers(pf, filename, center_cols, radius_col, unit='1'):
    """
    Return an iterator over EnzoSphere objects generated from the appropriate 
    columns in *filename*.  Optionally specify the *unit* radius is in.
    """
    sp_list = []
    for line in open(filename):
        if line.startswith("#"): continue
        vals = line.split()
        x,y,z = [float(vals[i]) for i in center_cols]
        r = float(vals[radius_col])
        yield pf.h.sphere([x,y,z], r/pf[unit])

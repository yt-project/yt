"""
A collection of helper functions, most generally for things
that SciPy doesn't have that I expected it to

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *
#import Interpolators as IT

class UnilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names):
        self.table = table
        x0, x1 = boundaries
        self.x_name = field_names
        self.x_bins = na.linspace(x0, x1, table.shape[0])

    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel()

        x_i = na.digitize(x_vals, self.x_bins) - 1
        if na.any((x_i == -1) | (x_i == len(self.x_bins)-1)):
            mylog.error("Sorry, but your values are outside" + \
                        " the table!  Dunno what to do, so dying.")
            mylog.error("Error was in: %s", data_object)
            raise ValueError

        x = (x_vals - self.x_bins[x_i]) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        xm = (self.x_bins[x_i+1] - x_vals) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        my_vals  = self.table[x_i  ] * (xm)
        my_vals += self.table[x_i+1] * (x )
        return my_vals.reshape(orig_shape)

class BilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names, truncate=False):
        self.table = table
        self.truncate = truncate
        x0, x1, y0, y1 = boundaries
        self.x_name, self.y_name = field_names
        self.x_bins = na.linspace(x0, x1, table.shape[0])
        self.y_bins = na.linspace(y0, y1, table.shape[1])

    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel()
        y_vals = data_object[self.y_name].ravel()

        x_i = na.digitize(x_vals, self.x_bins) - 1
        y_i = na.digitize(y_vals, self.y_bins) - 1
        if na.any((x_i == -1) | (x_i == len(self.x_bins)-1)) \
            or na.any((y_i == -1) | (y_i == len(self.y_bins)-1)):
            if not self.truncate:
                mylog.error("Sorry, but your values are outside" + \
                            " the table!  Dunno what to do, so dying.")
                mylog.error("Error was in: %s", data_object)
                raise ValueError
            else:
                x_i = na.minimum(na.maximum(x_i,0), len(self.x_bins)-2)
                y_i = na.minimum(na.maximum(y_i,0), len(self.y_bins)-2)

        my_vals = na.zeros(x_vals.shape, dtype='float')
        IT.BilinearlyInterpolate(self.table,
                                 x_vals, y_vals, self.x_bins, self.y_bins,
                                 x_i, y_i, my_vals)
        return my_vals

        x = (x_vals - self.x_bins[x_i]) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        y = (y_vals - self.y_bins[y_i]) / (self.y_bins[y_i+1] - self.y_bins[y_i])
        xm = (self.x_bins[x_i+1] - x_vals) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        ym = (self.y_bins[y_i+1] - y_vals) / (self.y_bins[y_i+1] - self.y_bins[y_i])
        my_vals  = self.table[x_i  ,y_i  ] * (xm*ym)
        my_vals += self.table[x_i+1,y_i  ] * (x *ym)
        my_vals += self.table[x_i  ,y_i+1] * (xm*y )
        my_vals += self.table[x_i+1,y_i+1] * (x *y )
        return my_vals.reshape(orig_shape)

class TrilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names, truncate = False):
        self.table = table
        self.truncate = truncate
        x0, x1, y0, y1, z0, z1 = boundaries
        self.x_name, self.y_name, self.z_name = field_names
        self.x_bins = na.linspace(x0, x1, table.shape[0])
        self.y_bins = na.linspace(y0, y1, table.shape[1])
        self.z_bins = na.linspace(z0, z1, table.shape[2])

    def __call__(self, data_object):
        orig_shape = data_object[self.x_name].shape
        x_vals = data_object[self.x_name].ravel()
        y_vals = data_object[self.y_name].ravel()
        z_vals = data_object[self.z_name].ravel()

        x_i = na.digitize(x_vals, self.x_bins) - 1
        y_i = na.digitize(y_vals, self.y_bins) - 1
        z_i = na.digitize(z_vals, self.z_bins) - 1
        if na.any((x_i == -1) | (x_i == len(self.x_bins)-1)) \
            or na.any((y_i == -1) | (y_i == len(self.y_bins)-1)) \
            or na.any((z_i == -1) | (z_i == len(self.z_bins)-1)):
            if not self.truncate:
                mylog.error("Sorry, but your values are outside" + \
                            " the table!  Dunno what to do, so dying.")
                mylog.error("Error was in: %s", data_object)
                raise ValueError
            else:
                x_i = na.minimum(na.maximum(x_i,0), len(self.x_bins)-2)
                y_i = na.minimum(na.maximum(y_i,0), len(self.y_bins)-2)
                z_i = na.minimum(na.maximum(z_i,0), len(self.z_bins)-2)

        my_vals = na.zeros(x_vals.shape, dtype='float')
        print self.table.dtype, x_vals.dtype, self.x_bins.dtype
        IT.TrilinearlyInterpolate(self.table,
                                 x_vals, y_vals, z_vals,
                                 self.x_bins, self.y_bins, self.z_bins,
                                 x_i, y_i, z_i, my_vals)
        return my_vals

        # Use notation from Paul Bourke's page on interpolation
        # http://local.wasp.uwa.edu.au/~pbourke/other/interpolation/
        x = (x_vals - self.x_bins[x_i]) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        y = (y_vals - self.y_bins[y_i]) / (self.y_bins[y_i+1] - self.y_bins[y_i])
        z = (z_vals - self.z_bins[z_i]) / (self.z_bins[z_i+1] - self.z_bins[z_i])
        xm = (self.x_bins[x_i+1] - x_vals) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        ym = (self.y_bins[y_i+1] - y_vals) / (self.y_bins[y_i+1] - self.y_bins[y_i])
        zm = (self.z_bins[z_i+1] - z_vals) / (self.z_bins[z_i+1] - self.z_bins[z_i])
        if na.any(na.isnan(self.table)):
            raise ValueError
        if na.any(na.isnan(x) | na.isnan(y) | na.isnan(z)):
            raise ValueError
        if na.any(na.isnan(xm) | na.isnan(ym) | na.isnan(zm)):
            raise ValueError
        my_vals  = self.table[x_i  ,y_i  ,z_i  ] * (xm*ym*zm)
        my_vals += self.table[x_i+1,y_i  ,z_i  ] * (x *ym*zm)
        my_vals += self.table[x_i  ,y_i+1,z_i  ] * (xm*y *zm)
        my_vals += self.table[x_i  ,y_i  ,z_i+1] * (xm*ym*z )
        my_vals += self.table[x_i+1,y_i  ,z_i+1] * (x *ym*z )
        my_vals += self.table[x_i  ,y_i+1,z_i+1] * (xm*y *z )
        my_vals += self.table[x_i+1,y_i+1,z_i  ] * (x *y *zm)
        my_vals += self.table[x_i+1,y_i+1,z_i+1] * (x *y *z )
        return my_vals.reshape(orig_shape)

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

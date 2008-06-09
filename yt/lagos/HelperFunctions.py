"""
A collection of helper functions, most generally for things
that SciPy doesn't have that I expected it to

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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
        my_vals = self.table[x_i  ] * (xm) \
                + self.table[x_i+1] * (x )
        return my_vals.reshape(orig_shape)

class BilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names):
        self.table = table
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
            mylog.error("Sorry, but your values are outside" + \
                        " the table!  Dunno what to do, so dying.")
            mylog.error("Error was in: %s", data_object)
            raise ValueError

        x = (x_vals - self.x_bins[x_i]) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        y = (y_vals - self.y_bins[y_i]) / (self.y_bins[y_i+1] - self.y_bins[y_i])
        xm = (self.x_bins[x_i+1] - x_vals) / (self.x_bins[x_i+1] - self.x_bins[x_i])
        ym = (self.y_bins[y_i+1] - y_vals) / (self.y_bins[y_i+1] - self.y_bins[y_i])
        my_vals = \
                  self.table[x_i  ,y_i  ] * (xm*ym) \
                + self.table[x_i+1,y_i  ] * (x *ym) \
                + self.table[x_i  ,y_i+1] * (xm*y ) \
                + self.table[x_i+1,y_i+1] * (x *y )
        return my_vals.reshape(orig_shape)

class TrilinearFieldInterpolator:
    def __init__(self, table, boundaries, field_names):
        self.table = table
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
            mylog.error("Sorry, but your values are outside" + \
                        " the table!  Dunno what to do, so dying.")
            mylog.error("Error was in: %s", data_object)
            raise ValueError

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
        my_vals = \
                  self.table[x_i  ,y_i  ,z_i  ] * (xm*ym*zm) \
                + self.table[x_i+1,y_i  ,z_i  ] * (x *ym*zm) \
                + self.table[x_i  ,y_i+1,z_i  ] * (xm*y *zm) \
                + self.table[x_i  ,y_i  ,z_i+1] * (xm*ym*z ) \
                + self.table[x_i+1,y_i  ,z_i+1] * (x *ym*z ) \
                + self.table[x_i  ,y_i+1,z_i+1] * (xm*y *z ) \
                + self.table[x_i+1,y_i+1,z_i  ] * (x *y *zm) \
                + self.table[x_i+1,y_i+1,z_i+1] * (x *y *z )
        return my_vals.reshape(orig_shape)

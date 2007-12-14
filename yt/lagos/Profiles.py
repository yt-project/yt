"""
Profile class, to deal with generating and obtaining profiles

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

# Note we do not inherit from EnzoData.
# We could, but I think we instead want to deal with the root datasource.
class BinnedProfile:
    def __init__(self, data_source, lazy_reader):
        self._data_source = data_source
        self._data = {}
        self._lazy_reader = lazy_reader

    def add_fields(self, fields, weight = "CellMassMsun", accumulation = False):
        # Note that the specification has to be the same for all of these
        fields = ensure_list(fields)
        if not self._lazy_reader:
            for field in fields:
                f, w, u = self._bin_field(self._data_source, field, weight,
                                          accumulation, self._args, check_cut = False)
                if weight:
                    f /= w
                self[field] = f
            self["UsedBins"] = u
        else:
            data = {}
            weight_data = {}
            for field in fields:
                data[field] = self._get_empty_field()
                weight_data[field] = self._get_empty_field()
            used = self._get_empty_field().astype('bool')
            for grid in self._data_source._grids:
                mylog.debug("Binner adding fields from grid %s", grid)
                args = self._get_bins(grid, check_cut=True)
                if not args:
                    continue
                for field in fields:
                    f, w, u = self._bin_field(grid, field, weight, accumulation,
                                              args=args, check_cut=True)
                    data[field] += f
                    weight_data[field] += w
                    used = (used | u)
                grid.clear_data()
            for field in fields:
                if weight:
                    data[field] /= weight_data[field]
                self[field] = data[field]
            self["UsedBins"] = used


    def __getitem__(self, key):
        # This raises a KeyError if it doesn't exist
        # This is because we explicitly want to
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value


# @todo: Fix accumulation with overriding
class BinnedProfile1D(BinnedProfile):
    def __init__(self, data_source, n_bins, lower_bound, upper_bound,
                 bin_field="Radius", log_space = True, lazy_reader=False):
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.bin_field = bin_field
        if log_space:
            func = na.logspace
            lower_bound, upper_bound = na.log10(lower_bound), na.log10(upper_bound)
        else:
            func = na.linspace
        self[bin_field] = func(lower_bound, upper_bound, n_bins)

        if not lazy_reader:
            self._args = self._get_bins(data_source)

    def _get_empty_field(self):
        return na.zeros(self[self.bin_field].size, dtype='float64')

    def _bin_field(self, source, field, weight, accumulation,
                   inv_bin_indices, check_cut=False):
        if check_cut:
            cm = self._data_source._get_point_indices(source)
            source_data = source[field][cm]
            if weight: weight_data = source[weight][cm]
        else:
            source_data = source[field]
            if weight: weight_data = source[weight]
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field() + 1.0
        used_field = na.ones(weight_field.shape, dtype='bool')
        for bin in inv_bin_indices.keys():
            temp_field = source_data[inv_bin_indices[bin]]
            if weight:
                weight_field[bin] = weight_data[inv_bin_indices[bin]].sum()
                temp_field *= weight_data[inv_bin_indices[bin]]
            binned_field[bin] = temp_field.sum()
        if accumulation: # Fix for laziness
            binned_field = na.add.accumulate(binned_field)
        return binned_field, weight_field

    def _get_bins(self, source, check_cut=False):
        if check_cut:
            cm = self._data_source._get_point_indices(source)
            source_data = source[self.bin_field][cm]
        else:
            source_data = source[self.bin_field]
        bin_order = na.argsort(source_data)
        bin_indices = na.searchsorted(self[self.bin_field],
                                      source_data)
        # Now we set up our inverse bin indices
        inv_bin_indices = {}
        for bin in range(self[self.bin_field].size):
            inv_bin_indices[bin] = na.where(bin_indices == bin)
        return inv_bin_indices


class BinnedProfile2D(BinnedProfile):
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 lazy_reader=False):

        self.total_cells = 0
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        if x_log:
            self[x_bin_field] = na.logspace(na.log10(x_lower_bound*0.99),
                                            na.log10(x_upper_bound*1.01),
                                            x_n_bins)
        else:
            self[x_bin_field] = na.linspace(
                x_lower_bound*0.99, x_upper_bound*1.01, x_n_bins)
        if y_log:
            self[y_bin_field] = na.logspace(na.log10(y_lower_bound*0.99),
                                            na.log10(y_upper_bound*1.01),
                                            y_n_bins)
        else:
            self[y_bin_field] = na.linspace(
                y_lower_bound*0.99, y_upper_bound*1.01, y_n_bins)
        if not lazy_reader:
            self._args = self._get_bins(data_source)
    def _get_empty_field(self):
        return na.zeros((self[self.x_bin_field].size,
                         self[self.y_bin_field].size), dtype='float64')

    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        #mylog.debug("Binning %s", field)
        if check_cut:
            pointI = self._data_source._get_point_indices(source)
            source_data = source[field][pointI].ravel()
            weight_data = na.ones(source_data.shape)
            if weight: weight_data = source[weight][pointI].ravel()
        else:
            source_data = source[field].ravel()
            weight_data = na.ones(source_data.shape)
            if weight: weight_data = source[weight].ravel()
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = na.ones(binned_field.shape,dtype='float64')
        used_field = self._get_empty_field()
        bin_indices_x = args[0].ravel()
        bin_indices_y = args[1].ravel()
        self.total_cells += source_data.size
        #mylog.debug("Binning %s times", source_data.size)
        nx = source_data.size
        try:
            weave.inline(_2d_profile_code, ['nx','bin_indices_x','bin_indices_y',
                                'weight_field','weight_data',
                                'binned_field','source_data', 'used_field'],
                         compiler='gcc', type_converters=converters.blitz,
                         auto_downcast = 0)
        except:
            mylog.debug("SciPy weaving failed; non-fatal, but slow.")
            for k in range(nx):
                i,j = bin_indices_x[k]-1, bin_indices_y[k]-1
                weight_field[i,j] += weight_data[k]
                binned_field[i,j] += source_data[k]*weight_data[k]
                used_field[i,j] = 1.0
        if accumulation: # Fix for laziness
            if not iterable(accumulation):
                raise SyntaxError("Accumulation needs to have length 2")
            if accumulation[0]:
                binned_field = na.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = na.add.accumulate(binned_field, axis=1)
        return binned_field, weight_field, used_field.astype('bool')

    def _get_bins(self, source, check_cut=False):
        if check_cut:
            cm = self._data_source._get_point_indices(source)
            source_data_x = source[self.x_bin_field][cm]
            source_data_y = source[self.y_bin_field][cm]
        else:
            source_data_x = source[self.x_bin_field]
            source_data_y = source[self.y_bin_field]
        if source_data_x.size == 0:
            return

        bin_indices_x = na.digitize(source_data_x.ravel(),
                                    self[self.x_bin_field])
        bin_indices_y = na.digitize(source_data_y.ravel(),
                                    self[self.y_bin_field])
        # Now we set up our inverse bin indices
        return (bin_indices_x, bin_indices_y)


_2d_profile_code = r"""
       int i,j;
       for(int n = 0; n < nx ; n++) {
         i = bin_indices_x(n)-1;
         j = bin_indices_y(n)-1;
         weight_field(i,j) += weight_data(n);
         binned_field(i,j) += source_data(n) * weight_data(n);
         used_field(i,j) = 1.0;
       }
       """
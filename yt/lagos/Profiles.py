"""
Profile classes, to deal with generating and obtaining profiles

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

def preserve_source_parameters(func):
    def save_state(*args, **kwargs):
        prof = args[0]
        source = args[1]
        if hasattr(source, 'field_parameters'):
            old_params = source.field_parameters
            source.field_parameters = prof._data_source.field_parameters
            tr = func(*args, **kwargs)
            source.field_parameters = old_params
        else:
            tr = func(*args, **kwargs)
        return tr
    return save_state

# Note we do not inherit from EnzoData.
# We could, but I think we instead want to deal with the root datasource.
class BinnedProfile:
    def __init__(self, data_source, lazy_reader):
        self._data_source = data_source
        self._data = {}
        self._lazy_reader = lazy_reader

    def _lazy_add_fields(self, fields, weight, accumulation):
        data = {}
        weight_data = {}
        for field in fields:
            data[field] = self._get_empty_field()
            weight_data[field] = self._get_empty_field()
        used = self._get_empty_field().astype('bool')
        pbar = get_pbar('Binning grids', len(self._data_source._grids))
        for gi,grid in enumerate(self._data_source._grids):
            pbar.update(gi)
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
        pbar.finish()
        ub = na.where(used)
        for field in fields:
            if weight:
                data[field][ub] /= weight_data[field][ub]
            self[field] = data[field]
        self["UsedBins"] = used

    def _unlazy_add_fields(self, fields, weight, accumulation):
        for field in fields:
            f, w, u = self._bin_field(self._data_source, field, weight,
                                      accumulation, self._args, check_cut = False)
            ub = na.where(u)
            if weight:
                f[ub] /= w[ub]
            self[field] = f
        self["UsedBins"] = u

    def add_fields(self, fields, weight = "CellMassMsun", accumulation = False):
        # Note that the specification has to be the same for all of these
        fields = ensure_list(fields)
        if self._lazy_reader:
            self._lazy_add_fields(fields, weight, accumulation)
        else:
            self._unlazy_add_fields(fields, weight, accumulation)

    def keys(self):
        return self._data.keys()

    def __getitem__(self, key):
        # This raises a KeyError if it doesn't exist
        # This is because we explicitly want to add all fields
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value


# @todo: Fix accumulation with overriding
class BinnedProfile1D(BinnedProfile):
    def __init__(self, data_source, n_bins, bin_field,
                 lower_bound, upper_bound,
                 log_space = True, lazy_reader=False):
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.bin_field = bin_field
        self._x_log = log_space
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

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        mi, inv_bin_indices = args
        if check_cut:
            cm = self._data_source._get_point_indices(source)
            source_data = source[field][cm].astype('float64')[mi]
            if weight: weight_data = source[weight][cm].astype('float64')[mi]
        else:
            source_data = source[field].astype('float64')[mi]
            if weight: weight_data = source[weight].astype('float64')[mi]
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        used_field = na.ones(weight_field.shape, dtype='bool')
        for bin in inv_bin_indices.keys():
            temp_field = source_data[inv_bin_indices[bin]]
            if weight:
                weight_field[bin] = weight_data[inv_bin_indices[bin]].sum()
                temp_field *= weight_data[inv_bin_indices[bin]]
            binned_field[bin] = temp_field.sum()
        if accumulation: # Fix for laziness
            binned_field = na.add.accumulate(binned_field)
        return binned_field, weight_field, na.ones(binned_field.shape,dtype='bool')

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        if check_cut:
            cm = self._data_source._get_point_indices(source)
            source_data = source[self.bin_field][cm]
        else:
            source_data = source[self.bin_field]
        if source_data.size == 0:
            return
        mi = na.where( (source_data > self[self.bin_field].min())
                     & (source_data < self[self.bin_field].max()))
        sd = source_data[mi]
        if sd.size == 0:
            return
        bin_indices = na.digitize(sd, self[self.bin_field])
        # Now we set up our inverse bin indices
        inv_bin_indices = {}
        for bin in range(self[self.bin_field].size):
            inv_bin_indices[bin] = na.where(bin_indices == bin)
        return (mi, inv_bin_indices)

    def write_out(self, filename, format="%0.16e"):
        fid = open(filename,"w")
        fields = [field for field in sorted(self._data.keys()) if field != "UsedBins"]
        fid.write("\t".join(["#"] + fields + ["\n"]))
        field_data = na.array([self._data[field] for field in fields])
        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

class BinnedProfile2D(BinnedProfile):
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 lazy_reader=False):

        self.total_cells = 0
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        self._x_log = x_log
        self._y_log = y_log
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
        if na.any(na.isnan(self[x_bin_field])) \
            or na.any(na.isnan(self[y_bin_field])):
            mylog.error("Your min/max values for x, y have given me a nan.")
            mylog.error("Usually this means you are asking for log, with a zero bound.")
            raise ValueError
        if not lazy_reader:
            self._args = self._get_bins(data_source)

    def _get_empty_field(self):
        return na.zeros((self[self.x_bin_field].size,
                         self[self.y_bin_field].size), dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        if check_cut:
            pointI = self._data_source._get_point_indices(source)
            source_data = source[field][pointI].ravel().astype('float64')
            weight_data = na.ones(source_data.shape).astype('float64')
            if weight: weight_data = source[weight][pointI].ravel().astype('float64')
        else:
            source_data = source[field].ravel().astype('float64')
            weight_data = na.ones(source_data.shape).astype('float64')
            if weight: weight_data = source[weight].ravel().astype('float64')
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        used_field = self._get_empty_field()
        mi = args[0]
        bin_indices_x = args[1].ravel().astype('int64')
        bin_indices_y = args[2].ravel().astype('int64')
        source_data = source_data[mi]
        weight_data = weight_data[mi]
        self.total_cells += bin_indices_x.size
        nx = bin_indices_x.size
        #mylog.debug("Binning %s / %s times", source_data.size, nx)
        PointCombine.Bin2DProfile(bin_indices_x, bin_indices_y, weight_data, source_data,
                     weight_field, binned_field, used_field)
        if accumulation: # Fix for laziness
            if not iterable(accumulation):
                raise SyntaxError("Accumulation needs to have length 2")
            if accumulation[0]:
                binned_field = na.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = na.add.accumulate(binned_field, axis=1)
        return binned_field, weight_field, used_field.astype('bool')

    @preserve_source_parameters
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
        mi = na.where( (source_data_x > self[self.x_bin_field].min())
                     & (source_data_x < self[self.x_bin_field].max())
                     & (source_data_y > self[self.y_bin_field].min())
                     & (source_data_y < self[self.y_bin_field].max()))
        sd_x = source_data_x[mi]
        sd_y = source_data_y[mi]
        if sd_x.size == 0 or sd_y.size == 0:
            return
        bin_indices_x = na.digitize(sd_x, self[self.x_bin_field])
        bin_indices_y = na.digitize(sd_y, self[self.y_bin_field])
        # Now we set up our inverse bin indices
        return (mi, bin_indices_x, bin_indices_y)

    def write_out(self, filename, format="%0.16e"):
        fid = open(filename,"w")
        fields = [field for field in sorted(self._data.keys()) if field != "UsedBins"]
        fid.write("\t".join(["#"] + [self.x_bin_field, self.y_bin_field]
                          + fields + ["\n"]))
        x,y = na.meshgrid(self._data[self.y_bin_field],
                          self._data[self.y_bin_field])
        field_data = [x.ravel(), y.ravel()]
        field_data += [self._data[field].ravel() for field in fields
                       if field not in [self.x_bin_field, self.y_bin_field]]
        field_data = na.array(field_data)
        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

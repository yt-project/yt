"""
Profile classes, to deal with generating and obtaining profiles

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

_field_mapping = {
    "total_mass": ("CellMassMsun", "ParticleMassMsun"),
    "hybrid_radius": ("RadiusCode", "ParticleRadiusCode"),
                 }

class EmptyProfileData(Exception):
    pass

def preserve_source_parameters(func):
    def save_state(*args, **kwargs):
        # Temporarily replace the 'field_parameters' for a
        # grid with the 'field_parameters' for the data source
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
class BinnedProfile(ParallelAnalysisInterface):
    def __init__(self, data_source, lazy_reader):
        self._data_source = data_source
        self.pf = data_source.pf
        self._data = {}
        self._pdata = {}
        self._lazy_reader = lazy_reader

    @property
    def hierarchy(self):
        return self.pf.hierarchy

    def _get_dependencies(self, fields):
        return ParallelAnalysisInterface._get_dependencies(
                    self, fields + self._get_bin_fields())

    def _initialize_parallel(self, fields):
        g_objs = [g for g in self._get_grid_objs()]
        self._preload(g_objs, self._get_dependencies(fields),
                      self._data_source.hierarchy.queue)

    def _lazy_add_fields(self, fields, weight, accumulation):
        self._ngrids = 0
        self.__data = {}         # final results will go here
        self.__weight_data = {}  # we need to track the weights as we go
        for field in fields:
            self.__data[field] = self._get_empty_field()
            self.__weight_data[field] = self._get_empty_field()
        self.__used = self._get_empty_field().astype('bool')
        #pbar = get_pbar('Binning grids', len(self._data_source._grids))
        for gi,grid in enumerate(self._get_grids(fields)):
            self._ngrids += 1
            #pbar.update(gi)
            args = self._get_bins(grid, check_cut=True)
            for field in fields:
                # We get back field values, weight values, used bins
                try:
                    f, w, u = self._bin_field(grid, field, weight, accumulation,
                                              args=args, check_cut=True)
                except EmptySourceData: 
                    # No bins returned for this grid, so forget it!
                    break
                self.__data[field] += f        # running total
                self.__weight_data[field] += w # running total
                self.__used = (self.__used | u)       # running 'or'
            grid.clear_data()
        # When the loop completes the parallel finalizer gets called
        #pbar.finish()
        ub = na.where(self.__used)
        for field in fields:
            if weight: # Now, at the end, we divide out.
                self.__data[field][ub] /= self.__weight_data[field][ub]
            self[field] = self.__data[field]
        self["UsedBins"] = self.__used
        del self.__data, self.__weight_data, self.__used

    def _finalize_parallel(self):
        for key in self.__data:
            self.__data[key] = self._mpi_allsum(self.__data[key])
        for key in self.__weight_data:
            self.__weight_data[key] = self._mpi_allsum(self.__weight_data[key])
        self.__used = self._mpi_allsum(self.__used)

    def _unlazy_add_fields(self, fields, weight, accumulation):
        for field in fields:
            f, w, u = self._bin_field(self._data_source, field, weight,
                                      accumulation, self._args, check_cut = False)
            if weight:
                f[u] /= w[u]
            self[field] = f
        self["myweight"] = w
        self["UsedBins"] = u

    def add_fields(self, fields, weight = "CellMassMsun", accumulation = False):
        """
        We accept a list of *fields* which will be binned if *weight* is not
        None and otherwise summed.  *accumulation* determines whether or not
        they will be accumulated from low to high along the appropriate axes.
        """
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

    def _get_field(self, source, this_field, check_cut):
        # This is where we will iterate to get all contributions to a field
        # which is how we will implement hybrid particle/cell fields
        # but...  we default to just the field.
        data = []
        for field in _field_mapping.get(this_field, (this_field,)):
            pointI = None
            if check_cut:
                if field in self.pf.field_info \
                    and self.pf.field_info[field].particle_type:
                    pointI = self._data_source._get_particle_indices(source)
                else:
                    pointI = self._data_source._get_point_indices(source)
            data.append(source[field][pointI].ravel().astype('float64'))
        return na.concatenate(data, axis=0)

    def _fix_pickle(self):
        if isinstance(self._data_source, tuple):
            self._data_source = self._data_source[1]

# @todo: Fix accumulation with overriding
class BinnedProfile1D(BinnedProfile):
    def __init__(self, data_source, n_bins, bin_field,
                 lower_bound, upper_bound,
                 log_space = True, lazy_reader=False,
                 left_collect = False):
        """
        A 'Profile' produces either a weighted (or unweighted) average or a
        straight sum of a field in a bin defined by another field.  In the case
        of a weighted average, we have: p_i = sum( w_i * v_i ) / sum(w_i)

        We accept a *data_source*, which will be binned into *n_bins* by the
        field *bin_field* between the *lower_bound* and the *upper_bound*.
        These bins may or may not be equally divided in *log_space*, and the
        *lazy_reader* flag controls whether we use a memory conservative
        approach.
        """
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.bin_field = bin_field
        self._x_log = log_space
        self.left_collect = left_collect
        # Get our bins
        if log_space:
            func = na.logspace
            lower_bound, upper_bound = na.log10(lower_bound), na.log10(upper_bound)
        else:
            func = na.linspace
        # These are the bin *edges*
        self._bins = func(lower_bound, upper_bound, n_bins + 1)

        # These are the bin *left edges*
        self[bin_field] = self._bins[:-1]

        # If we are not being memory-conservative, grab all the bins
        # and the inverse indices right now.
        if not lazy_reader:
            self._args = self._get_bins(data_source)

    def _get_empty_field(self):
        return na.zeros(self[self.bin_field].size, dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        mi, inv_bin_indices = args # Args has the indices to use as input
        # check_cut is set if source != self._data_source
        # (i.e., lazy_reader)
        source_data = self._get_field(source, field, check_cut)[mi]
        if weight: weight_data = self._get_field(source, weight, check_cut)[mi]
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        used_field = na.ones(weight_field.shape, dtype='bool')
        # Now we perform the actual binning
        for bin in inv_bin_indices.keys():
            # temp_field is *all* the points from source that go into this bin
            temp_field = source_data[inv_bin_indices[bin]]
            if weight:
                # now w_i * v_i and store sum(w_i)
                weight_field[bin] = weight_data[inv_bin_indices[bin]].sum()
                temp_field *= weight_data[inv_bin_indices[bin]]
            binned_field[bin] = temp_field.sum()
        # Fix for laziness, because at the *end* we will be
        # summing up all of the histograms and dividing by the
        # weights.  Accumulation likely doesn't work with weighted
        # average fields.
        if accumulation: 
            binned_field = na.add.accumulate(binned_field)
        return binned_field, weight_field, na.ones(binned_field.shape,dtype='bool')

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data = self._get_field(source, self.bin_field, check_cut)
        if source_data.size == 0: # Nothing for us here.
            raise EmptyProfileData()
        # Truncate at boundaries.
        if self.left_collect:
            mi = na.where(source_data < self._bins.max())
        else:
            mi = na.where( (source_data > self._bins.min())
                         & (source_data < self._bins.max()))
        sd = source_data[mi]
        if sd.size == 0:
            raise EmptyProfileData()
        # Stick the bins into our fixed bins, set at initialization
        bin_indices = na.digitize(sd, self._bins) - 1
        if self.left_collect: bin_indices = na.maximum(0, bin_indices)
        # Now we set up our inverse bin indices
        inv_bin_indices = {}
        for bin in range(self[self.bin_field].size):
            # Which fall into our bin?
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

    def _get_bin_fields(self):
        return [self.bin_field]

class BinnedProfile2D(BinnedProfile):
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 lazy_reader=False, left_collect=False):
        """
        A 'Profile' produces either a weighted (or unweighted) average or a
        straight sum of a field in a bin defined by two other fields.  In the case
        of a weighted average, we have: p_i = sum( w_i * v_i ) / sum(w_i)

        We accept a *data_source*, which will be binned into *x_n_bins* by the
        field *x_bin_field* between the *x_lower_bound* and the *x_upper_bound*
        and then again binned into *y_n_bins* by the field *y_bin_field*
        between the *y_lower_bound* and the *y_upper_bound*.  These bins may or
        may not be equally divided in log-space as specified by *x_log* and *y_log*,
        and the *lazy_reader* flag controls whether we use a memory
        conservative approach.
        """
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        self._x_log = x_log
        self._y_log = y_log
        self.left_collect = left_collect

        func = {True:na.logspace, False:na.linspace}[x_log]
        bounds = fix_bounds(x_upper_bound, x_lower_bound, x_log)
        self._x_bins = func(bounds[0], bounds[1], x_n_bins + 1)
        self[x_bin_field] = self._x_bins[:-1]

        func = {True:na.logspace, False:na.linspace}[y_log]
        bounds = fix_bounds(y_upper_bound, y_lower_bound, y_log)
        self._y_bins = func(bounds[0], bounds[1], y_n_bins + 1)
        self[y_bin_field] = self._y_bins[:-1]


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
        source_data = self._get_field(source, field, check_cut)
        if weight: weight_data = self._get_field(source, weight, check_cut)
        else: weight_data = na.ones(source_data.shape, dtype='float64')
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        used_field = self._get_empty_field()
        mi = args[0]
        bin_indices_x = args[1].ravel().astype('int64')
        bin_indices_y = args[2].ravel().astype('int64')
        source_data = source_data[mi]
        weight_data = weight_data[mi]
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
        source_data_x = self._get_field(source, self.x_bin_field, check_cut)
        source_data_y = self._get_field(source, self.y_bin_field, check_cut)
        if source_data_x.size == 0:
            raise EmptyProfileData()
        if self.left_collect:
            mi = na.where( (source_data_x < self._x_bins.max())
                         & (source_data_y < self._y_bins.max()))
        else:
            mi = na.where( (source_data_x > self._x_bins.min())
                         & (source_data_x < self._x_bins.max())
                         & (source_data_y > self._y_bins.min())
                         & (source_data_y < self._y_bins.max()))
        sd_x = source_data_x[mi]
        sd_y = source_data_y[mi]
        if sd_x.size == 0 or sd_y.size == 0:
            raise EmptyProfileData()
        bin_indices_x = na.digitize(sd_x, self._x_bins) - 1
        bin_indices_y = na.digitize(sd_y, self._y_bins) - 1
        if self.left_collect:
            bin_indices_x = na.maximum(bin_indices_x, 0)
            bin_indices_y = na.maximum(bin_indices_y, 0)
        # Now we set up our inverse bin indices
        return (mi, bin_indices_x, bin_indices_y)

    def write_out(self, filename, format="%0.16e"):
        """
        Write out the values of x,y,v in ascii to *filename* for every field in
        the profile.  Optionally a *format* can be specified.
        """
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

    def _get_bin_fields(self):
        return [self.x_bin_field, self.y_bin_field]

def fix_bounds(upper, lower, logit):
    if logit: return na.log10(upper), na.log10(lower)
    return upper, lower

class BinnedProfile3D(BinnedProfile):
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 z_n_bins, z_bin_field, z_lower_bound, z_upper_bound, z_log,
                 lazy_reader=False):
        BinnedProfile.__init__(self, data_source, lazy_reader)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        self.z_bin_field = z_bin_field
        self._x_log = x_log
        self._y_log = y_log
        self._z_log = z_log

        func = {True:na.logspace, False:na.linspace}[x_log]
        bounds = fix_bounds(x_upper_bound, x_lower_bound, x_log)
        self._x_bins = func(bounds[0], bounds[1], x_n_bins + 1)
        self[x_bin_field] = self._x_bins[:-1]

        func = {True:na.logspace, False:na.linspace}[y_log]
        bounds = fix_bounds(y_upper_bound, y_lower_bound, y_log)
        self._y_bins = func(bounds[0], bounds[1], y_n_bins + 1)
        self[y_bin_field] = self._y_bins[:-1]

        func = {True:na.logspace, False:na.linspace}[z_log]
        bounds = fix_bounds(z_upper_bound, z_lower_bound, z_log)
        self._z_bins = func(bounds[0], bounds[1], z_n_bins + 1)
        self[z_bin_field] = self._z_bins[:-1]

        if na.any(na.isnan(self[x_bin_field])) \
            or na.any(na.isnan(self[y_bin_field])) \
            or na.any(na.isnan(self[z_bin_field])):
            mylog.error("Your min/max values for x, y or z have given me a nan.")
            mylog.error("Usually this means you are asking for log, with a zero bound.")
            raise ValueError
        if not lazy_reader:
            self._args = self._get_bins(data_source)

    def _get_empty_field(self):
        return na.zeros((self[self.x_bin_field].size,
                         self[self.y_bin_field].size,
                         self[self.z_bin_field].size), dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        source_data = self._get_field(source, field, check_cut)
        weight_data = na.ones(source_data.shape).astype('float64')
        if weight: weight_data = self._get_field(source, weight, check_cut)
        else: weight_data = na.ones(source_data.shape).astype('float64')
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        used_field = self._get_empty_field()
        mi = args[0]
        bin_indices_x = args[1].ravel().astype('int64')
        bin_indices_y = args[2].ravel().astype('int64')
        bin_indices_z = args[3].ravel().astype('int64')
        source_data = source_data[mi]
        weight_data = weight_data[mi]
        PointCombine.Bin3DProfile(
                     bin_indices_x, bin_indices_y, bin_indices_z,
                     weight_data, source_data,
                     weight_field, binned_field, used_field)
        if accumulation: # Fix for laziness
            if not iterable(accumulation):
                raise SyntaxError("Accumulation needs to have length 2")
            if accumulation[0]:
                binned_field = na.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = na.add.accumulate(binned_field, axis=1)
            if accumulation[2]:
                binned_field = na.add.accumulate(binned_field, axis=2)
        return binned_field, weight_field, used_field.astype('bool')

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data_x = self._get_field(source, self.x_bin_field, check_cut)
        source_data_y = self._get_field(source, self.y_bin_field, check_cut)
        source_data_y = self._get_field(source, self.z_bin_field, check_cut)
        if source_data_x.size == 0:
            raise EmptyProfileData()
        mi = ( (source_data_x > self._x_bins.min())
             & (source_data_x < self._x_bins.max())
             & (source_data_y > self._y_bins.min())
             & (source_data_y < self._y_bins.max())
             & (source_data_z > self._z_bins.min())
             & (source_data_z < self._z_bins.max()))
        sd_x = source_data_x[mi]
        sd_y = source_data_y[mi]
        sd_z = source_data_z[mi]
        if sd_x.size == 0 or sd_y.size == 0 or sd_z.size == 0:
            raise EmptyProfileData()
        bin_indices_x = na.digitize(sd_x, self._x_bins)
        bin_indices_y = na.digitize(sd_y, self._y_bins)
        bin_indices_z = na.digitize(sd_z, self._z_bins)
        # Now we set up our inverse bin indices
        return (mi, bin_indices_x, bin_indices_y, bin_indices_z)

    def write_out(self, filename, format="%0.16e"):
        pass # Will eventually dump HDF5

    def _get_bin_fields(self):
        return [self.x_bin_field, self.y_bin_field, self.z_bin_field]

    def store_profile(self, name, force=False):
        """
        By identifying the profile with a fixed, user-input *name* we can
        store it in the serialized data section of the hierarchy file.  *force*
        governs whether or not an existing profile with that name will be
        overwritten.
        """
        # First we get our data in order
        order = []
        set_attr = {'x_bin_field':self.x_bin_field,
                    'y_bin_field':self.y_bin_field,
                    'z_bin_field':self.z_bin_field,
                    'x_bin_values':self[self.x_bin_field],
                    'y_bin_values':self[self.y_bin_field],
                    'z_bin_values':self[self.z_bin_field],
                    '_x_log':self._x_log,
                    '_y_log':self._y_log,
                    '_z_log':self._z_log,
                    'shape': (self[self.x_bin_field].size,
                              self[self.y_bin_field].size,
                              self[self.z_bin_field].size),
                    'field_order':order }
        values = []
        for field in self._data:
            if field in set_attr.values(): continue
            order.append(field)
            values.append(self[field].ravel())
        values = na.array(values).transpose()
        self._data_source.hierarchy.save_data(values, "/Profiles", name,
                                              set_attr, force=force)

class StoredBinnedProfile3D(BinnedProfile3D):
    def __init__(self, pf, name):
        """
        Given a *pf* parameterfile and the *name* of a stored profile, retrieve
        it into a read-only data structure.
        """
        self._data = {}
        prof_arr = pf.h.get_data("/Profiles", name)
        if prof_arr is None: raise KeyError("No such array")
        for ax in 'xyz':
            for base in ['%s_bin_field', '_%s_log']:
                setattr(self, base % ax, prof_arr.getAttr(base % ax))
        for ax in 'xyz':
            fn = getattr(self, '%s_bin_field' % ax)
            self._data[fn] = prof_arr.getAttr('%s_bin_values' % ax)
        shape = prof_arr.getAttr('shape')
        for fn, fd in zip(prof_arr.getAttr('field_order'),
                          prof_arr.read().transpose()):
            self._data[fn] = fd.reshape(shape)

    def add_fields(self, *args, **kwargs):
        raise RuntimeError("Sorry, you can't add to a stored profile.")

"""
Profile classes, to deal with generating and obtaining profiles



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np

from yt.funcs import *

from yt.units.yt_array import uconcatenate, array_like_field
from yt.units.unit_object import Unit
from yt.data_objects.data_containers import YTFieldData
from yt.utilities.lib.misc_utilities import \
    bin_profile1d, bin_profile2d, bin_profile3d, \
    new_bin_profile1d, new_bin_profile2d, \
    new_bin_profile3d
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_objects
from yt.utilities.exceptions import YTEmptyProfileData

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
    def __init__(self, data_source):
        ParallelAnalysisInterface.__init__(self)
        self._data_source = data_source
        self.ds = data_source.ds
        self.field_data = YTFieldData()

    @property
    def index(self):
        return self.ds.index

    def _get_dependencies(self, fields):
        return ParallelAnalysisInterface._get_dependencies(
                    self, fields + self._get_bin_fields())

    def add_fields(self, fields, weight = "cell_mass", accumulation = False, fractional=False):
        """
        We accept a list of *fields* which will be binned if *weight* is not
        None and otherwise summed.  *accumulation* determines whether or not
        they will be accumulated from low to high along the appropriate axes.
        """
        # Note that the specification has to be the same for all of these
        fields = ensure_list(fields)
        data = {}         # final results will go here
        weight_data = {}  # we need to track the weights as we go
        std_data = {}
        for field in fields:
            data[field] = self._get_empty_field()
            weight_data[field] = self._get_empty_field()
            std_data[field] = self._get_empty_field()
        used = self._get_empty_field().astype('bool')
        chunk_fields = fields[:]
        if weight is not None: chunk_fields += [weight]
        #pbar = get_pbar('Binning grids', len(self._data_source._grids))
        for ds in self._data_source.chunks(chunk_fields, chunking_style = "io"):
            try:
                args = self._get_bins(ds, check_cut=True)
            except YTEmptyProfileData:
                # No bins returned for this grid, so forget it!
                continue
            for field in fields:
                # We get back field values, weight values, used bins
                f, w, q, u = self._bin_field(ds, field, weight, accumulation,
                                          args=args, check_cut=True)
                data[field] += f        # running total
                weight_data[field] += w # running total
                used |= u       # running 'or'
                std_data[field][u] += w[u] * (q[u]/w[u] + \
                    (f[u]/w[u] -
                     data[field][u]/weight_data[field][u])**2) # running total
        for key in data:
            data[key] = self.comm.mpi_allreduce(data[key], op='sum')
        for key in weight_data:
            weight_data[key] = self.comm.mpi_allreduce(weight_data[key], op='sum')
        used = self.comm.mpi_allreduce(used, op='sum')
        # When the loop completes the parallel finalizer gets called
        #pbar.finish()
        ub = np.where(used)
        for field in fields:
            if weight: # Now, at the end, we divide out.
                data[field][ub] /= weight_data[field][ub]
                std_data[field][ub] /= weight_data[field][ub]
            self[field] = data[field]
            self["%s_std" % field] = np.sqrt(std_data[field])
        self["UsedBins"] = used

        if fractional:
            for field in fields:
                self.field_data[field] /= self.field_data[field].sum()

    def keys(self):
        return self.field_data.keys()

    def __getitem__(self, key):
        # This raises a KeyError if it doesn't exist
        # This is because we explicitly want to add all fields
        return self.field_data[key]

    def __setitem__(self, key, value):
        self.field_data[key] = value

    def _get_field(self, source, field, check_cut):
        # This is where we will iterate to get all contributions to a field
        # which is how we will implement hybrid particle/cell fields
        # but...  we default to just the field.
        data = []
        data.append(source[field].astype('float64'))
        return uconcatenate(data, axis=0)

    def _fix_pickle(self):
        if isinstance(self._data_source, tuple):
            self._data_source = self._data_source[1]

# @todo: Fix accumulation with overriding
class BinnedProfile1D(BinnedProfile):
    """
    A 'Profile' produces either a weighted (or unweighted) average or a
    straight sum of a field in a bin defined by another field.  In the case
    of a weighted average, we have: p_i = sum( w_i * v_i ) / sum(w_i)

    We accept a *data_source*, which will be binned into *n_bins*
    by the field *bin_field* between the *lower_bound* and the
    *upper_bound*.  These bins may or may not be equally divided
    in *log_space*, and the *lazy_reader* flag controls whether we
    use a memory conservative approach. If *end_collect* is True,
    take all values outside the given bounds and store them in the
    0 and *n_bins*-1 values.
    """
    def __init__(self, data_source, n_bins, bin_field,
                 lower_bound, upper_bound,
                 log_space = True,
                 end_collect=False):
        BinnedProfile.__init__(self, data_source)
        self.bin_field = bin_field
        self._x_log = log_space
        self.end_collect = end_collect
        self.n_bins = n_bins

        # Get our bins
        if log_space:
            if lower_bound <= 0.0 or upper_bound <= 0.0:
                raise YTIllDefinedBounds(lower_bound, upper_bound)
            func = np.logspace
            lower_bound, upper_bound = np.log10(lower_bound), np.log10(upper_bound)
        else:
            func = np.linspace

        # These are the bin *edges*
        self._bins = func(lower_bound, upper_bound, n_bins + 1)

        # These are the bin *left edges*.  These are the x-axis values
        # we plot in the PlotCollection
        self[bin_field] = self._bins

        # If we are not being memory-conservative, grab all the bins
        # and the inverse indices right now.

    def _get_empty_field(self):
        return np.zeros(self[self.bin_field].size, dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        mi, inv_bin_indices = args # Args has the indices to use as input
        # check_cut is set if source != self._data_source
        source_data = self._get_field(source, field, check_cut)
        if weight: weight_data = self._get_field(source, weight, check_cut)
        else: weight_data = np.ones(source_data.shape, dtype='float64')
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        m_field = self._get_empty_field()
        q_field = self._get_empty_field()
        used_field = self._get_empty_field()
        mi = args[0]
        bin_indices_x = args[1].ravel().astype('int64')
        source_data = source_data[mi]
        weight_data = weight_data[mi]
        bin_profile1d(bin_indices_x, weight_data, source_data,
                      weight_field, binned_field,
                      m_field, q_field, used_field)
        # Fix for laziness, because at the *end* we will be
        # summing up all of the histograms and dividing by the
        # weights.  Accumulation likely doesn't work with weighted
        # average fields.
        if accumulation:
            binned_field = np.add.accumulate(binned_field)
        return binned_field, weight_field, q_field, \
            used_field.astype("bool")

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data = self._get_field(source, self.bin_field, check_cut)
        if source_data.size == 0: # Nothing for us here.
            raise YTEmptyProfileData()
        # Truncate at boundaries.
        if self.end_collect:
            mi = np.ones_like(source_data).astype('bool')
        else:
            mi = ((source_data > self._bins.min())
               &  (source_data < self._bins.max()))
        sd = source_data[mi]
        if sd.size == 0:
            raise YTEmptyProfileData()
        # Stick the bins into our fixed bins, set at initialization
        bin_indices = np.digitize(sd, self._bins)
        if self.end_collect: #limit the range of values to 0 and n_bins-1
            bin_indices = np.clip(bin_indices, 0, self.n_bins - 1)
        else: #throw away outside values
            bin_indices -= 1

        return (mi, bin_indices)

    def choose_bins(self, bin_style):
        # Depending on the bin_style, choose from bin edges 0...N either:
        # both: 0...N, left: 0...N-1, right: 1...N
        # center: N bins that are the average (both in linear or log
        # space) of each pair of left/right edges
        x = self.field_data[self.bin_field]
        if bin_style is 'both': pass
        elif bin_style is 'left': x = x[:-1]
        elif bin_style is 'right': x = x[1:]
        elif bin_style is 'center':
            if self._x_log: x=np.log10(x)
            x = 0.5*(x[:-1] + x[1:])
            if self._x_log: x=10**x
        else:
            mylog.error('Did not recognize bin_style')
            raise ValueError
        return x

    def write_out(self, filename, format="%0.16e", bin_style='left'):
        '''
        Write out data in ascii file, using *format* and
        *bin_style* (left, right, center, both).
        '''
        fid = open(filename,"w")
        fields = [field for field in sorted(self.field_data.keys()) if field != "UsedBins"]
        fields.remove(self.bin_field)
        fid.write("\t".join(["#"] + [self.bin_field] + fields + ["\n"]))

        field_data = np.array(self.choose_bins(bin_style))
        if bin_style is 'both':
            field_data = np.append([field_data], np.array([self.field_data[field] for field in fields]), axis=0)
        else:
            field_data = np.append([field_data], np.array([self.field_data[field][:-1] for field in fields]), axis=0)

        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

    def write_out_h5(self, filename, group_prefix=None, bin_style='left'):
        """
        Write out data in an hdf5 file *filename*.  Each profile is
        put into a group, named by the axis fields.  Optionally a
        *group_prefix* can be prepended to the group name.  If the
        group already exists, it will delete and replace.  However,
        due to hdf5 functionality, in only unlinks the data, so an
        h5repack may be necessary to conserve space.  Axes values are
        saved in group attributes.  Bins will be saved based on
        *bin_style* (left, right, center, both).
        """
        fid = h5py.File(filename)
        fields = [field for field in sorted(self.field_data.keys()) if (field != "UsedBins" and field != self.bin_field)]
        if group_prefix is None:
            name = "%s-1d" % (self.bin_field)
        else:
            name = "%s-%s-1d" % (group_prefix, self.bin_field)

        if name in fid:
            mylog.info("Profile file is getting larger since you are attempting to overwrite a profile. You may want to repack")
            del fid[name]
        group = fid.create_group(name)
        group.attrs["x-axis-%s" % self.bin_field] = self.choose_bins(bin_style)
        for field in fields:
            dset = group.create_dataset("%s" % field, data=self.field_data[field][:-1])
        fid.close()

    def _get_bin_fields(self):
        return [self.bin_field]

class BinnedProfile2D(BinnedProfile):
    """
    A 'Profile' produces either a weighted (or unweighted) average
    or a straight sum of a field in a bin defined by two other
    fields.  In the case of a weighted average, we have: p_i =
    sum( w_i * v_i ) / sum(w_i)

    We accept a *data_source*, which will be binned into
    *x_n_bins* by the field *x_bin_field* between the
    *x_lower_bound* and the *x_upper_bound* and then again binned
    into *y_n_bins* by the field *y_bin_field* between the
    *y_lower_bound* and the *y_upper_bound*.  These bins may or
    may not be equally divided in log-space as specified by
    *x_log* and *y_log*, and the *lazy_reader* flag controls
    whether we use a memory conservative approach. If
    *end_collect* is True, take all values outside the given
    bounds and store them in the 0 and *n_bins*-1 values.
    """
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 end_collect=False):
        BinnedProfile.__init__(self, data_source)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        self._x_log = x_log
        self._y_log = y_log
        self.end_collect = end_collect
        self.x_n_bins = x_n_bins
        self.y_n_bins = y_n_bins

        func = {True:np.logspace, False:np.linspace}[x_log]
        bounds = fix_bounds(x_lower_bound, x_upper_bound, x_log)
        self._x_bins = func(bounds[0], bounds[1], x_n_bins + 1)
        self[x_bin_field] = self._x_bins

        func = {True:np.logspace, False:np.linspace}[y_log]
        bounds = fix_bounds(y_lower_bound, y_upper_bound, y_log)
        self._y_bins = func(bounds[0], bounds[1], y_n_bins + 1)
        self[y_bin_field] = self._y_bins

        if np.any(np.isnan(self[x_bin_field])) \
            or np.any(np.isnan(self[y_bin_field])):
            mylog.error("Your min/max values for x, y have given me a nan.")
            mylog.error("Usually this means you are asking for log, with a zero bound.")
            raise ValueError

    def _get_empty_field(self):
        return np.zeros((self[self.x_bin_field].size,
                         self[self.y_bin_field].size), dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        source_data = self._get_field(source, field, check_cut)
        if weight: weight_data = self._get_field(source, weight, check_cut)
        else: weight_data = np.ones(source_data.shape, dtype='float64')
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        m_field = self._get_empty_field()
        q_field = self._get_empty_field()
        used_field = self._get_empty_field()
        mi = args[0]
        bin_indices_x = args[1].ravel().astype('int64')
        bin_indices_y = args[2].ravel().astype('int64')
        source_data = source_data[mi]
        weight_data = weight_data[mi]
        nx = bin_indices_x.size
        #mylog.debug("Binning %s / %s times", source_data.size, nx)
        bin_profile2d(bin_indices_x, bin_indices_y, weight_data, source_data,
                      weight_field, binned_field, m_field, q_field, used_field)
        if accumulation: # Fix for laziness
            if not iterable(accumulation):
                raise SyntaxError("Accumulation needs to have length 2")
            if accumulation[0]:
                binned_field = np.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = np.add.accumulate(binned_field, axis=1)
        return binned_field, weight_field, q_field, \
            used_field.astype("bool")

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data_x = self._get_field(source, self.x_bin_field, check_cut)
        source_data_y = self._get_field(source, self.y_bin_field, check_cut)
        if source_data_x.size == 0:
            raise YTEmptyProfileData()

        if self.end_collect:
            mi = np.arange(source_data_x.size)
        else:
            mi = np.where( (source_data_x > self._x_bins.min())
                           & (source_data_x < self._x_bins.max())
                           & (source_data_y > self._y_bins.min())
                           & (source_data_y < self._y_bins.max()))
        sd_x = source_data_x[mi]
        sd_y = source_data_y[mi]
        if sd_x.size == 0 or sd_y.size == 0:
            raise YTEmptyProfileData()

        bin_indices_x = np.digitize(sd_x, self._x_bins) - 1
        bin_indices_y = np.digitize(sd_y, self._y_bins) - 1
        if self.end_collect:
            bin_indices_x = np.minimum(np.maximum(1, bin_indices_x), self.x_n_bins) - 1
            bin_indices_y = np.minimum(np.maximum(1, bin_indices_y), self.y_n_bins) - 1

        # Now we set up our inverse bin indices
        return (mi, bin_indices_x, bin_indices_y)

    def choose_bins(self, bin_style):
        # Depending on the bin_style, choose from bin edges 0...N either:
        # both: 0...N, left: 0...N-1, right: 1...N
        # center: N bins that are the average (both in linear or log
        # space) of each pair of left/right edges

        x = self.field_data[self.x_bin_field]
        y = self.field_data[self.y_bin_field]
        if bin_style is 'both':
            pass
        elif bin_style is 'left':
            x = x[:-1]
            y = y[:-1]
        elif bin_style is 'right':
            x = x[1:]
            y = y[1:]
        elif bin_style is 'center':
            if self._x_log: x=np.log10(x)
            if self._y_log: y=np.log10(y)
            x = 0.5*(x[:-1] + x[1:])
            y = 0.5*(y[:-1] + y[1:])
            if self._x_log: x=10**x
            if self._y_log: y=10**y
        else:
            mylog.error('Did not recognize bin_style')
            raise ValueError

        return x,y

    def write_out(self, filename, format="%0.16e", bin_style='left'):
        """
        Write out the values of x,y,v in ascii to *filename* for every
        field in the profile.  Optionally a *format* can be specified.
        Bins will be saved based on *bin_style* (left, right, center,
        both).
        """
        fid = open(filename,"w")
        fields = [field for field in sorted(self.field_data.keys()) if field != "UsedBins"]
        fid.write("\t".join(["#"] + [self.x_bin_field, self.y_bin_field]
                          + fields + ["\n"]))
        x,y = self.choose_bins(bin_style)
        x,y = np.meshgrid(x,y)
        field_data = [x.ravel(), y.ravel()]
        if bin_style is not 'both':
            field_data += [self.field_data[field][:-1,:-1].ravel() for field in fields
                           if field not in [self.x_bin_field, self.y_bin_field]]
        else:
            field_data += [self.field_data[field].ravel() for field in fields
                           if field not in [self.x_bin_field, self.y_bin_field]]

        field_data = np.array(field_data)
        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

    def write_out_h5(self, filename, group_prefix=None, bin_style='left'):
        """
        Write out data in an hdf5 file.  Each profile is put into a
        group, named by the axis fields.  Optionally a group_prefix
        can be prepended to the group name.  If the group already
        exists, it will delete and replace.  However, due to hdf5
        functionality, in only unlinks the data, so an h5repack may be
        necessary to conserve space.  Axes values are saved in group
        attributes. Bins will be saved based on *bin_style* (left,
        right, center, both).
        """
        fid = h5py.File(filename)
        fields = [field for field in sorted(self.field_data.keys()) if (field != "UsedBins" and field != self.x_bin_field and field != self.y_bin_field)]
        if group_prefix is None:
            name = "%s-%s-2d" % (self.y_bin_field, self.x_bin_field)
        else:
            name = "%s-%s-%s-2d" % (group_prefix, self.y_bin_field, self.x_bin_field)
        if name in fid:
            mylog.info("Profile file is getting larger since you are attempting to overwrite a profile. You may want to repack")
            del fid[name]
        group = fid.create_group(name)

        xbins, ybins = self.choose_bins(bin_style)
        group.attrs["x-axis-%s" % self.x_bin_field] = xbins
        group.attrs["y-axis-%s" % self.y_bin_field] = ybins
        for field in fields:
            dset = group.create_dataset("%s" % field, data=self.field_data[field][:-1,:-1])
        fid.close()

    def _get_bin_fields(self):
        return [self.x_bin_field, self.y_bin_field]

def fix_bounds(upper, lower, logit):
    if logit:
        if lower <= 0.0 or upper <= 0.0:
            raise YTIllDefinedBounds(lower, upper)
        return np.log10(upper), np.log10(lower)
    return upper, lower

class BinnedProfile3D(BinnedProfile):
    """
    A 'Profile' produces either a weighted (or unweighted) average
    or a straight sum of a field in a bin defined by two other
    fields.  In the case of a weighted average, we have: p_i =
    sum( w_i * v_i ) / sum(w_i)

    We accept a *data_source*, which will be binned into
    *(x,y,z)_n_bins* by the field *(x,y,z)_bin_field* between the
    *(x,y,z)_lower_bound* and the *(x,y,z)_upper_bound*.  These bins may or
    may not be equally divided in log-space as specified by *(x,y,z)_log*.
    If *end_collect* is True, take all values outside the given bounds and
    store them in the 0 and *n_bins*-1 values.
    """
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 z_n_bins, z_bin_field, z_lower_bound, z_upper_bound, z_log,
                 end_collect=False):
        BinnedProfile.__init__(self, data_source)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        self.z_bin_field = z_bin_field
        self._x_log = x_log
        self._y_log = y_log
        self._z_log = z_log
        self.end_collect = end_collect
        self.x_n_bins = x_n_bins
        self.y_n_bins = y_n_bins
        self.z_n_bins = z_n_bins

        func = {True:np.logspace, False:np.linspace}[x_log]
        bounds = fix_bounds(x_lower_bound, x_upper_bound, x_log)
        self._x_bins = func(bounds[0], bounds[1], x_n_bins + 1)
        self[x_bin_field] = self._x_bins

        func = {True:np.logspace, False:np.linspace}[y_log]
        bounds = fix_bounds(y_lower_bound, y_upper_bound, y_log)
        self._y_bins = func(bounds[0], bounds[1], y_n_bins + 1)
        self[y_bin_field] = self._y_bins

        func = {True:np.logspace, False:np.linspace}[z_log]
        bounds = fix_bounds(z_lower_bound, z_upper_bound, z_log)
        self._z_bins = func(bounds[0], bounds[1], z_n_bins + 1)
        self[z_bin_field] = self._z_bins

        if np.any(np.isnan(self[x_bin_field])) \
            or np.any(np.isnan(self[y_bin_field])) \
            or np.any(np.isnan(self[z_bin_field])):
            mylog.error("Your min/max values for x, y or z have given me a nan.")
            mylog.error("Usually this means you are asking for log, with a zero bound.")
            raise ValueError

    def _get_empty_field(self):
        return np.zeros((self[self.x_bin_field].size,
                         self[self.y_bin_field].size,
                         self[self.z_bin_field].size), dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        source_data = self._get_field(source, field, check_cut)
        weight_data = np.ones(source_data.shape).astype('float64')
        if weight: weight_data = self._get_field(source, weight, check_cut)
        else: weight_data = np.ones(source_data.shape).astype('float64')
        self.total_stuff = source_data.sum()
        binned_field = self._get_empty_field()
        weight_field = self._get_empty_field()
        m_field = self._get_empty_field()
        q_field = self._get_empty_field()
        used_field = self._get_empty_field()
        mi = args[0]
        bin_indices_x = args[1].ravel().astype('int64')
        bin_indices_y = args[2].ravel().astype('int64')
        bin_indices_z = args[3].ravel().astype('int64')
        source_data = source_data[mi]
        weight_data = weight_data[mi]
        bin_profile3d(bin_indices_x, bin_indices_y, bin_indices_z,
                      weight_data, source_data, weight_field, binned_field,
                      m_field, q_field, used_field)
        if accumulation: # Fix for laziness
            if not iterable(accumulation):
                raise SyntaxError("Accumulation needs to have length 2")
            if accumulation[0]:
                binned_field = np.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = np.add.accumulate(binned_field, axis=1)
            if accumulation[2]:
                binned_field = np.add.accumulate(binned_field, axis=2)
        return binned_field, weight_field, q_field, \
            used_field.astype("bool")

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data_x = self._get_field(source, self.x_bin_field, check_cut)
        source_data_y = self._get_field(source, self.y_bin_field, check_cut)
        source_data_z = self._get_field(source, self.z_bin_field, check_cut)
        if source_data_x.size == 0:
            raise YTEmptyProfileData()
        if self.end_collect:
            mi = np.arange(source_data_x.size)
        else:
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
            raise YTEmptyProfileData()

        bin_indices_x = np.digitize(sd_x, self._x_bins) - 1
        bin_indices_y = np.digitize(sd_y, self._y_bins) - 1
        bin_indices_z = np.digitize(sd_z, self._z_bins) - 1
        if self.end_collect:
            bin_indices_x = np.minimum(np.maximum(1, bin_indices_x), self.x_n_bins) - 1
            bin_indices_y = np.minimum(np.maximum(1, bin_indices_y), self.y_n_bins) - 1
            bin_indices_z = np.minimum(np.maximum(1, bin_indices_z), self.z_n_bins) - 1

        # Now we set up our inverse bin indices
        return (mi, bin_indices_x, bin_indices_y, bin_indices_z)

    def choose_bins(self, bin_style):
        # Depending on the bin_style, choose from bin edges 0...N either:
        # both: 0...N, left: 0...N-1, right: 1...N
        # center: N bins that are the average (both in linear or log
        # space) of each pair of left/right edges

        x = self.field_data[self.x_bin_field]
        y = self.field_data[self.y_bin_field]
        z = self.field_data[self.z_bin_field]
        if bin_style is 'both':
            pass
        elif bin_style is 'left':
            x = x[:-1]
            y = y[:-1]
            z = z[:-1]
        elif bin_style is 'right':
            x = x[1:]
            y = y[1:]
            z = z[1:]
        elif bin_style is 'center':
            if self._x_log: x=np.log10(x)
            if self._y_log: y=np.log10(y)
            if self._z_log: z=np.log10(z)
            x = 0.5*(x[:-1] + x[1:])
            y = 0.5*(y[:-1] + y[1:])
            z = 0.5*(z[:-1] + z[1:])
            if self._x_log: x=10**x
            if self._y_log: y=10**y
            if self._z_log: y=10**z
        else:
            mylog.error('Did not recognize bin_style')
            raise ValueError

        return x,y,z

    def write_out(self, filename, format="%0.16e"):
        pass # Will eventually dump HDF5

    def write_out_h5(self, filename, group_prefix=None, bin_style='left'):
        """
        Write out data in an hdf5 file.  Each profile is put into a
        group, named by the axis fields.  Optionally a group_prefix
        can be prepended to the group name.  If the group already
        exists, it will delete and replace.  However, due to hdf5
        functionality, in only unlinks the data, so an h5repack may be
        necessary to conserve space.  Axes values are saved in group
        attributes.
        """
        fid = h5py.File(filename)
        fields = [field for field in sorted(self.field_data.keys())
                  if (field != "UsedBins" and field != self.x_bin_field and field != self.y_bin_field and field != self.z_bin_field)]
        if group_prefix is None:
            name = "%s-%s-%s-3d" % (self.z_bin_field, self.y_bin_field, self.x_bin_field)
        else:
            name = "%s-%s-%s-%s-3d" % (group_prefix,self.z_bin_field, self.y_bin_field, self.x_bin_field)

        if name in fid:
            mylog.info("Profile file is getting larger since you are attempting to overwrite a profile. You may want to repack")
            del fid[name]
        group = fid.create_group(name)

        xbins, ybins, zbins= self.choose_bins(bin_style)
        group.attrs["x-axis-%s" % self.x_bin_field] = xbins
        group.attrs["y-axis-%s" % self.y_bin_field] = ybins
        group.attrs["z-axis-%s" % self.z_bin_field] = zbins

        for field in fields:
            dset = group.create_dataset("%s" % field, data=self.field_data[field][:-1,:-1,:-1])
        fid.close()


    def _get_bin_fields(self):
        return [self.x_bin_field, self.y_bin_field, self.z_bin_field]

    def store_profile(self, name, force=False):
        """
        By identifying the profile with a fixed, user-input *name* we can
        store it in the serialized data section of the index file.  *force*
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
        for field in self.field_data:
            if field in set_attr.values(): continue
            order.append(field)
            values.append(self[field].ravel())
        values = np.array(values).transpose()
        self._data_source.index.save_data(values, "/Profiles", name,
                                              set_attr, force=force)

class ProfileFieldAccumulator(object):
    def __init__(self, n_fields, size):
        shape = size + (n_fields,)
        self.values = np.zeros(shape, dtype="float64")
        self.mvalues = np.zeros(shape, dtype="float64")
        self.qvalues = np.zeros(shape, dtype="float64")
        self.used = np.zeros(size, dtype='bool')
        self.weight_values = np.zeros(size, dtype="float64")

class ProfileND(ParallelAnalysisInterface):
    """The profile object class"""
    def __init__(self, data_source, weight_field = None):
        self.data_source = data_source
        self.ds = data_source.ds
        self.field_map = {}
        self.field_data = YTFieldData()
        if weight_field is not None:
            self.variance = YTFieldData()
            weight_field = self.data_source._determine_fields(weight_field)[0]
        self.weight_field = weight_field
        self.field_units = {}
        ParallelAnalysisInterface.__init__(self, comm=data_source.comm)

    def add_fields(self, fields):
        """Add fields to profile

        Parameters
        ----------
        fields : list of field names
            A list of fields to create profile histograms for
        
        """
        fields = self.data_source._determine_fields(fields)
        temp_storage = ProfileFieldAccumulator(len(fields), self.size)
        cfields = fields + list(self.bin_fields)
        citer = self.data_source.chunks(cfields, "io")
        for chunk in parallel_objects(citer):
            self._bin_chunk(chunk, fields, temp_storage)
        self._finalize_storage(fields, temp_storage)

    def set_field_unit(self, field, new_unit):
        """Sets a new unit for the requested field

        Parameters
        ----------
        field : string or field tuple
           The name of the field that is to be changed.

        new_unit : string or Unit object
           The name of the new unit.
        """
        if field in self.field_units:
            self.field_units[field] = \
                Unit(new_unit, registry=self.ds.unit_registry)
        else:
            fd = self.field_map[field]
            if fd in self.field_units:
                self.field_units[fd] = \
                    Unit(new_unit, registry=self.ds.unit_registry)
            else:
                raise KeyError("%s not in profile!" % (field))

    def _finalize_storage(self, fields, temp_storage):
        # We use our main comm here
        # This also will fill _field_data

        for i, field in enumerate(fields):
            # q values are returned as q * weight but we want just q
            temp_storage.qvalues[..., i][temp_storage.used] /= \
              temp_storage.weight_values[temp_storage.used]

        # get the profile data from all procs
        all_store = {self.comm.rank: temp_storage}
        all_store = self.comm.par_combine_object(all_store,
                                                 "join", datatype="dict")

        all_val = np.zeros_like(temp_storage.values)
        all_mean = np.zeros_like(temp_storage.mvalues)
        all_var = np.zeros_like(temp_storage.qvalues)
        all_weight = np.zeros_like(temp_storage.weight_values)
        all_used = np.zeros_like(temp_storage.used, dtype="bool")

        # Combine the weighted mean and variance from each processor.
        # For two samples with total weight, mean, and variance 
        # given by w, m, and s, their combined mean and variance are:
        # m12 = (m1 * w1 + m2 * w2) / (w1 + w2)
        # s12 = (m1 * (s1**2 + (m1 - m12)**2) + 
        #        m2 * (s2**2 + (m2 - m12)**2)) / (w1 + w2)
        # Here, the mvalues are m and the qvalues are s**2.
        for p in sorted(all_store.keys()):
            all_used += all_store[p].used
            old_mean = all_mean.copy()
            old_weight = all_weight.copy()
            all_weight[all_store[p].used] += \
              all_store[p].weight_values[all_store[p].used]
            for i, field in enumerate(fields):
                all_val[..., i][all_store[p].used] += \
                  all_store[p].values[..., i][all_store[p].used]

                all_mean[..., i][all_store[p].used] = \
                  (all_mean[..., i] * old_weight +
                   all_store[p].mvalues[..., i] *
                   all_store[p].weight_values)[all_store[p].used] / \
                   all_weight[all_store[p].used]

                all_var[..., i][all_store[p].used] = \
                  (old_weight * (all_var[..., i] +
                                 (old_mean[..., i] - all_mean[..., i])**2) +
                   all_store[p].weight_values *
                   (all_store[p].qvalues[..., i] + 
                    (all_store[p].mvalues[..., i] -
                     all_mean[..., i])**2))[all_store[p].used] / \
                    all_weight[all_store[p].used]

        all_var = np.sqrt(all_var)
        del all_store
        self.used = all_used
        blank = ~all_used

        self.weight = all_weight
        self.weight[blank] = 0.0

        for i, field in enumerate(fields):
            if self.weight_field is None:
                self.field_data[field] = \
                  array_like_field(self.data_source, 
                                   all_val[...,i], field)
            else:
                self.field_data[field] = \
                  array_like_field(self.data_source, 
                                   all_mean[...,i], field)
                self.variance[field] = \
                  array_like_field(self.data_source,
                                   all_var[...,i], field)
                self.variance[field][blank] = 0.0
            self.field_data[field][blank] = 0.0
            self.field_units[field] = self.field_data[field].units
            if isinstance(field, tuple):
                self.field_map[field[1]] = field
            else:
                self.field_map[field] = field

    def _bin_chunk(self, chunk, fields, storage):
        raise NotImplementedError

    def _filter(self, bin_fields):
        # cut_points is set to be everything initially, but
        # we also want to apply a filtering based on min/max
        filter = np.ones(bin_fields[0].shape, dtype='bool')
        for (mi, ma), data in zip(self.bounds, bin_fields):
            filter &= (data > mi)
            filter &= (data < ma)
        return filter, [data[filter] for data in bin_fields]

    def _get_data(self, chunk, fields):
        # We are using chunks now, which will manage the field parameters and
        # the like.
        bin_fields = [chunk[bf] for bf in self.bin_fields]
        # We want to make sure that our fields are within the bounds of the
        # binning
        filter, bin_fields = self._filter(bin_fields)
        if not np.any(filter): return None
        arr = np.zeros((bin_fields[0].size, len(fields)), dtype="float64")
        for i, field in enumerate(fields):
            units = chunk.ds.field_info[field].units
            arr[:,i] = chunk[field][filter].in_units(units)
        if self.weight_field is not None:
            units = chunk.ds.field_info[self.weight_field].units
            weight_data = chunk[self.weight_field].in_units(units)
        else:
            weight_data = np.ones(filter.size, dtype="float64")
        weight_data = weight_data[filter]
        # So that we can pass these into
        return arr, weight_data, bin_fields

    def __getitem__(self, field):
        fname = self.field_map.get(field, None)
        if fname is None and isinstance(field, tuple):
            fname = self.field_map.get(field[1], None)
        if fname is None:
            raise KeyError(field)
        else:
            if getattr(self, 'fractional', False):
                return self.field_data[fname]
            else:
                return self.field_data[fname].in_units(self.field_units[fname])

    def items(self):
        return [(k,self[k]) for k in self.field_data.keys()]

    def __iter__(self):
        return sorted(self.items())

    def _get_bins(self, mi, ma, n, take_log):
        if take_log:
            return np.logspace(np.log10(mi), np.log10(ma), n+1)
        else:
            return np.linspace(mi, ma, n+1)

class Profile1D(ProfileND):
    """An object that represents a 1D profile.

    Parameters
    ----------

    data_source : AMD3DData object
        The data object to be profiled
    x_field : string field name
        The field to profile as a function of
    x_n : integer
        The number of bins along the x direction.
    x_min : float
        The minimum value of the x profile field.
    x_max : float
        The maximum value of the x profile field.
    x_log : boolean
        Controls whether or not the bins for the x field are evenly
        spaced in linear (False) or log (True) space.
    weight_field : string field name
        The field to weight the profiled fields by.

    """
    def __init__(self, data_source, x_field, x_n, x_min, x_max, x_log,
                 weight_field = None):
        super(Profile1D, self).__init__(data_source, weight_field)
        self.x_field = x_field
        self.x_log = x_log
        self.x_bins = array_like_field(data_source,
                                       self._get_bins(x_min, x_max, x_n, x_log),
                                       self.x_field)
        self.size = (self.x_bins.size - 1,)
        self.bin_fields = (self.x_field,)
        self.x = 0.5*(self.x_bins[1:]+self.x_bins[:-1])

    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None: return
        fdata, wdata, (bf_x,) = rv
        bin_ind = np.digitize(bf_x, self.x_bins) - 1
        new_bin_profile1d(bin_ind, wdata, fdata,
                      storage.weight_values, storage.values,
                      storage.mvalues, storage.qvalues,
                      storage.used)

        # We've binned it!

    def set_x_unit(self, new_unit):
        """Sets a new unit for the x field

        parameters
        ----------
        new_unit : string or Unit object
           The name of the new unit.
        """
        self.x_bins.convert_to_units(new_unit)
        self.x = 0.5*(self.x_bins[1:]+self.x_bins[:-1])

    @property
    def bounds(self):
        return ((self.x_bins[0], self.x_bins[-1]),)

class Profile2D(ProfileND):
    """An object that represents a 2D profile.

    Parameters
    ----------

    data_source : AMD3DData object
        The data object to be profiled
    x_field : string field name
        The field to profile as a function of along the x axis.
    x_n : integer
        The number of bins along the x direction.
    x_min : float
        The minimum value of the x profile field.
    x_max : float
        The maximum value of hte x profile field.
    x_log : boolean
        Controls whether or not the bins for the x field are evenly
        spaced in linear (False) or log (True) space.
    y_field : string field name
        The field to profile as a function of along the y axis
    y_n : integer
        The number of bins along the y direction.
    y_min : float
        The minimum value of the y profile field.
    y_max : float
        The maximum value of hte y profile field.
    y_log : boolean
        Controls whether or not the bins for the y field are evenly
        spaced in linear (False) or log (True) space.
    weight_field : string field name
        The field to weight the profiled fields by.

    """
    def __init__(self, data_source,
                 x_field, x_n, x_min, x_max, x_log,
                 y_field, y_n, y_min, y_max, y_log,
                 weight_field = None):
        super(Profile2D, self).__init__(data_source, weight_field)
        # X
        self.x_field = x_field
        self.x_log = x_log
        self.x_bins = array_like_field(data_source,
                                       self._get_bins(x_min, x_max, x_n, x_log),
                                       self.x_field)
        # Y
        self.y_field = y_field
        self.y_log = y_log
        self.y_bins = array_like_field(data_source,
                                       self._get_bins(y_min, y_max, y_n, y_log),
                                       self.y_field)

        self.size = (self.x_bins.size - 1, self.y_bins.size - 1)

        self.bin_fields = (self.x_field, self.y_field)
        self.x = 0.5*(self.x_bins[1:]+self.x_bins[:-1])
        self.y = 0.5*(self.y_bins[1:]+self.y_bins[:-1])

    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None: return
        fdata, wdata, (bf_x, bf_y) = rv
        bin_ind_x = np.digitize(bf_x, self.x_bins) - 1
        bin_ind_y = np.digitize(bf_y, self.y_bins) - 1
        new_bin_profile2d(bin_ind_x, bin_ind_y, wdata, fdata,
                      storage.weight_values, storage.values,
                      storage.mvalues, storage.qvalues,
                      storage.used)
        # We've binned it!

    def set_x_unit(self, new_unit):
        """Sets a new unit for the x field

        parameters
        ----------
        new_unit : string or Unit object
           The name of the new unit.
        """
        self.x_bins.convert_to_units(new_unit)
        self.x = 0.5*(self.x_bins[1:]+self.x_bins[:-1])

    def set_y_unit(self, new_unit):
        """Sets a new unit for the y field

        parameters
        ----------
        new_unit : string or Unit object
           The name of the new unit.
        """
        self.y_bins.convert_to_units(new_unit)
        self.y = 0.5*(self.y_bins[1:]+self.y_bins[:-1])

    @property
    def bounds(self):
        return ((self.x_bins[0], self.x_bins[-1]),
                (self.y_bins[0], self.y_bins[-1]))

class Profile3D(ProfileND):
    """An object that represents a 2D profile.

    Parameters
    ----------

    data_source : AMD3DData object
        The data object to be profiled
    x_field : string field name
        The field to profile as a function of along the x axis.
    x_n : integer
        The number of bins along the x direction.
    x_min : float
        The minimum value of the x profile field.
    x_max : float
        The maximum value of hte x profile field.
    x_log : boolean
        Controls whether or not the bins for the x field are evenly
        spaced in linear (False) or log (True) space.
    y_field : string field name
        The field to profile as a function of along the y axis
    y_n : integer
        The number of bins along the y direction.
    y_min : float
        The minimum value of the y profile field.
    y_max : float
        The maximum value of hte y profile field.
    y_log : boolean
        Controls whether or not the bins for the y field are evenly
        spaced in linear (False) or log (True) space.
    z_field : string field name
        The field to profile as a function of along the z axis
    z_n : integer
        The number of bins along the z direction.
    z_min : float
        The minimum value of the z profile field.
    z_max : float
        The maximum value of hte z profile field.
    z_log : boolean
        Controls whether or not the bins for the z field are evenly
        spaced in linear (False) or log (True) space.
    weight_field : string field name
        The field to weight the profiled fields by.

    """
    def __init__(self, data_source,
                 x_field, x_n, x_min, x_max, x_log,
                 y_field, y_n, y_min, y_max, y_log,
                 z_field, z_n, z_min, z_max, z_log,
                 weight_field = None):
        super(Profile3D, self).__init__(data_source, weight_field)
        # X
        self.x_field = x_field
        self.x_log = x_log
        self.x_bins = array_like_field(data_source,
                                       self._get_bins(x_min, x_max, x_n, x_log),
                                       self.x_field)
        # Y
        self.y_field = y_field
        self.y_log = y_log
        self.y_bins = array_like_field(data_source,
                                       self._get_bins(y_min, y_max, y_n, y_log),
                                       self.y_field)
        # Z
        self.z_field = z_field
        self.z_log = z_log
        self.z_bins = array_like_field(data_source,
                                       self._get_bins(z_min, z_max, z_n, z_log),
                                       self.z_field)

        self.size = (self.x_bins.size - 1,
                     self.y_bins.size - 1,
                     self.z_bins.size - 1)

        self.bin_fields = (self.x_field, self.y_field, self.z_field)
        self.x = 0.5*(self.x_bins[1:]+self.x_bins[:-1])
        self.y = 0.5*(self.y_bins[1:]+self.y_bins[:-1])
        self.z = 0.5*(self.z_bins[1:]+self.z_bins[:-1])

    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None: return
        fdata, wdata, (bf_x, bf_y, bf_z) = rv
        bin_ind_x = np.digitize(bf_x, self.x_bins) - 1
        bin_ind_y = np.digitize(bf_y, self.y_bins) - 1
        bin_ind_z = np.digitize(bf_z, self.z_bins) - 1
        new_bin_profile3d(bin_ind_x, bin_ind_y, bin_ind_z, wdata, fdata,
                      storage.weight_values, storage.values,
                      storage.mvalues, storage.qvalues,
                      storage.used)
        # We've binned it!

    @property
    def bounds(self):
        return ((self.x_bins[0], self.x_bins[-1]),
                (self.y_bins[0], self.y_bins[-1]),
                (self.z_bins[0], self.z_bins[-1]))

    def set_x_unit(self, new_unit):
        """Sets a new unit for the x field

        parameters
        ----------
        new_unit : string or Unit object
           The name of the new unit.
        """
        self.x_bins.convert_to_units(new_unit)
        self.x = 0.5*(self.x_bins[1:]+self.x_bins[:-1])

    def set_y_unit(self, new_unit):
        """Sets a new unit for the y field

        parameters
        ----------
        new_unit : string or Unit object
           The name of the new unit.
        """
        self.y_bins.convert_to_units(new_unit)
        self.y = 0.5*(self.y_bins[1:]+self.y_bins[:-1])

    def set_z_unit(self, new_unit):
        """Sets a new unit for the z field

        parameters
        ----------
        new_unit : string or Unit object
           The name of the new unit.
        """
        self.z_bins.convert_to_units(new_unit)
        self.z = 0.5*(self.z_bins[1:]+self.z_bins[:-1])


def sanitize_field_tuple_keys(input_dict, data_source):
    if input_dict is not None:
        dummy = {}
        for item in input_dict:
            dummy[data_source._determine_fields(item)[0]] = input_dict[item]
        return dummy
    else:
        return input_dict

def create_profile(data_source, bin_fields, fields, n_bins=64,
                   extrema=None, logs=None, units=None,
                   weight_field="cell_mass",
                   accumulation=False, fractional=False):
    r"""
    Create a 1, 2, or 3D profile object.

    The dimensionality of the profile object is chosen by the number of
    fields given in the bin_fields argument.

    Parameters
    ----------
    data_source : YTSelectionContainer Object
        The data object to be profiled.
    bin_fields : list of strings
        List of the binning fields for profiling.
    fields : list of strings
        The fields to be profiled.
    n_bins : int or list of ints
        The number of bins in each dimension.  If None, 64 bins for
        each bin are used for each bin field.
        Default: 64.
    extrema : dict of min, max tuples
        Minimum and maximum values of the bin_fields for the profiles.
        The keys correspond to the field names. Defaults to the extrema
        of the bin_fields of the dataset. If a units dict is provided, extrema
        are understood to be in the units specified in the dictionary.
    logs : dict of boolean values
        Whether or not to log the bin_fields for the profiles.
        The keys correspond to the field names. Defaults to the take_log
        attribute of the field.
    units : dict of strings
        The units of the fields in the profiles, including the bin_fields.
    weight_field : str or tuple field identifier
        The weight field for computing weighted average for the profile
        values.  If None, the profile values are sums of the data in
        each bin.
    accumulation : bool or list of bools
        If True, the profile values for a bin n are the cumulative sum of
        all the values from bin 0 to n.  If -True, the sum is reversed so
        that the value for bin n is the cumulative sum from bin N (total bins)
        to n.  If the profile is 2D or 3D, a list of values can be given to
        control the summation in each dimension independently.
        Default: False.
    fractional : If True the profile values are divided by the sum of all
        the profile data such that the profile represents a probability
        distribution function.

    Examples
    --------

    Create a 1d profile.  Access bin field from profile.x and field
    data from profile[<field_name>].

    >>> ds = load("DD0046/DD0046")
    >>> ad = ds.h.all_data()
    >>> profile = create_profile(ad, [("gas", "density")], 
    ...                              [("gas", "temperature"),
    ...                               ("gas", "velocity_x")])
    >>> print profile.x
    >>> print profile["gas", "temperature"]

    """
    bin_fields = data_source._determine_fields(bin_fields)
    fields = ensure_list(fields)
    if len(bin_fields) == 1:
        cls = Profile1D
    elif len(bin_fields) == 2:
        cls = Profile2D
    elif len(bin_fields) == 3:
        cls = Profile3D
    else:
        raise NotImplementedError
    bin_fields = data_source._determine_fields(bin_fields)
    fields = data_source._determine_fields(fields)
    units = sanitize_field_tuple_keys(units, data_source)
    extrema = sanitize_field_tuple_keys(extrema, data_source)
    logs = sanitize_field_tuple_keys(logs, data_source)
    if weight_field is not None:
        weight_field, = data_source._determine_fields([weight_field])
    if not iterable(n_bins):
        n_bins = [n_bins] * len(bin_fields)
    if not iterable(accumulation):
        accumulation = [accumulation] * len(bin_fields)
    if logs is None:
        logs = {}
    logs_list = []
    for bin_field in bin_fields:
        if bin_field in logs:
            logs_list.append(logs[bin_field])
        else:
            logs_list.append(data_source.ds.field_info[bin_field].take_log)
    logs = logs_list
    if extrema is None:
        ex = [data_source.quantities["Extrema"](f, non_zero=l)
              for f, l in zip(bin_fields, logs)]
    else:
        ex = []
        for bin_field in bin_fields:
            bf_units = data_source.ds.field_info[bin_field].units
            try:
                field_ex = list(extrema[bin_field[-1]])
            except KeyError:
                field_ex = list(extrema[bin_field])
            if units is not None and bin_field in units:
                if isinstance(field_ex[0], tuple):
                    field_ex = [data_source.ds.quan(*f) for f in field_ex]
                fe = data_source.ds.arr(field_ex, units[bin_field])
                fe.convert_to_units(bf_units)
                field_ex = [fe[0].v, fe[1].v]
            if iterable(field_ex[0]):
                field_ex[0] = data_source.ds.quan(field_ex[0][0], field_ex[0][1])
                field_ex[0] = field_ex[0].in_units(bf_units)
            if iterable(field_ex[1]):
                field_ex[1] = data_source.ds.quan(field_ex[1][0], field_ex[1][1])
                field_ex[1] = field_ex[1].in_units(bf_units)
            ex.append(field_ex)
    args = [data_source]
    for f, n, (mi, ma), l in zip(bin_fields, n_bins, ex, logs):
        args += [f, n, mi, ma, l]
    obj = cls(*args, weight_field = weight_field)
    setattr(obj, "accumulation", accumulation)
    setattr(obj, "fractional", fractional)
    if fields is not None:
        obj.add_fields([field for field in fields])
    for field in fields:
        if fractional:
            obj.field_data[field] /= obj.field_data[field].sum()
        for axis, acc in enumerate(accumulation):
            if not acc: continue
            temp = obj.field_data[field]
            temp = np.rollaxis(temp, axis)
            if weight_field is not None:
                temp_weight = obj.weight
                temp_weight = np.rollaxis(temp_weight, axis)
            if acc < 0:
                temp = temp[::-1]
                if weight_field is not None:
                    temp_weight = temp_weight[::-1]
            if weight_field is None:
                temp = temp.cumsum(axis=0)
            else:
                temp = (temp * temp_weight).cumsum(axis=0) / \
                  temp_weight.cumsum(axis=0)
            if acc < 0:
                temp = temp[::-1]
                if weight_field is not None:
                    temp_weight = temp_weight[::-1]
            temp = np.rollaxis(temp, axis)
            obj.field_data[field] = temp
            if weight_field is not None:
                temp_weight = np.rollaxis(temp_weight, axis)
                obj.weight = temp_weight
    if units is not None:
        for field, unit in units.iteritems():
            field = data_source._determine_fields(field)[0]
            if field == obj.x_field:
                obj.set_x_unit(unit)
            elif field == getattr(obj, "y_field", None):
                obj.set_y_unit(unit)
            elif field == getattr(obj, "z_field", None):
                obj.set_z_unit(unit)
            else:
                obj.set_field_unit(field, unit)
    return obj


"""
Profile classes, to deal with generating and obtaining profiles

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Samuel Skillman <samskillman@gmail.com>
Affiliation: CASA, University of Colorado at Boulder
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

import h5py
import numpy as na

from yt.funcs import *

from yt.data_objects.data_containers import YTFieldData
from yt.utilities.lib import bin_profile1d, bin_profile2d, bin_profile3d
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface

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
    def __init__(self, data_source):
        ParallelAnalysisInterface.__init__(self)
        self._data_source = data_source
        self.pf = data_source.pf
        self.field_data = YTFieldData()
        self._pdata = {}

    @property
    def hierarchy(self):
        return self.pf.hierarchy

    def _get_dependencies(self, fields):
        return ParallelAnalysisInterface._get_dependencies(
                    self, fields + self._get_bin_fields())

    def add_fields(self, fields, weight = "CellMassMsun", accumulation = False, fractional=False):
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
            except EmptyProfileData:
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
        ub = na.where(used)
        for field in fields:
            if weight: # Now, at the end, we divide out.
                data[field][ub] /= weight_data[field][ub]
                std_data[field][ub] /= weight_data[field][ub]
            self[field] = data[field]
            #self["%s_std" % field] = na.sqrt(std_data[field])
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

    def _get_field(self, source, this_field, check_cut):
        # This is where we will iterate to get all contributions to a field
        # which is how we will implement hybrid particle/cell fields
        # but...  we default to just the field.
        data = []
        for field in _field_mapping.get(this_field, (this_field,)):
            data.append(source[field].astype('float64'))
        return na.concatenate(data, axis=0)

    def _fix_pickle(self):
        if isinstance(self._data_source, tuple):
            self._data_source = self._data_source[1]

# @todo: Fix accumulation with overriding
class BinnedProfile1D(BinnedProfile):
    def __init__(self, data_source, n_bins, bin_field,
                 lower_bound, upper_bound,
                 log_space = True, 
                 end_collect=False):
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
        BinnedProfile.__init__(self, data_source)
        self.bin_field = bin_field
        self._x_log = log_space
        self.end_collect = end_collect
        self.n_bins = n_bins

        # Get our bins
        if log_space:
            func = na.logspace
            lower_bound, upper_bound = na.log10(lower_bound), na.log10(upper_bound)
        else:
            func = na.linspace

        # These are the bin *edges*
        self._bins = func(lower_bound, upper_bound, n_bins + 1)

        # These are the bin *left edges*.  These are the x-axis values
        # we plot in the PlotCollection
        self[bin_field] = self._bins

        # If we are not being memory-conservative, grab all the bins
        # and the inverse indices right now.

    def _get_empty_field(self):
        return na.zeros(self[self.bin_field].size, dtype='float64')

    @preserve_source_parameters
    def _bin_field(self, source, field, weight, accumulation,
                   args, check_cut=False):
        mi, inv_bin_indices = args # Args has the indices to use as input
        # check_cut is set if source != self._data_source
        # (i.e., lazy_reader)
        source_data = self._get_field(source, field, check_cut)
        if weight: weight_data = self._get_field(source, weight, check_cut)
        else: weight_data = na.ones(source_data.shape, dtype='float64')
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
            binned_field = na.add.accumulate(binned_field)
        return binned_field, weight_field, q_field, \
            used_field.astype("bool")

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data = self._get_field(source, self.bin_field, check_cut)
        if source_data.size == 0: # Nothing for us here.
            raise EmptyProfileData()
        # Truncate at boundaries.
        if self.end_collect:
            mi = na.ones_like(source_data).astype('bool')
        else:
            mi = ((source_data > self._bins.min())
               &  (source_data < self._bins.max()))
        sd = source_data[mi]
        if sd.size == 0:
            raise EmptyProfileData()
        # Stick the bins into our fixed bins, set at initialization
        bin_indices = na.digitize(sd, self._bins)
        if self.end_collect: #limit the range of values to 0 and n_bins-1
            bin_indices = na.clip(bin_indices, 0, self.n_bins - 1)
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
            if self._x_log: x=na.log10(x)
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

        field_data = na.array(self.choose_bins(bin_style)) 
        if bin_style is 'both':
            field_data = na.append([field_data], na.array([self.field_data[field] for field in fields]), axis=0)
        else: 
            field_data = na.append([field_data], na.array([self.field_data[field][:-1] for field in fields]), axis=0)
        
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
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 end_collect=False):
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
        BinnedProfile.__init__(self, data_source)
        self.x_bin_field = x_bin_field
        self.y_bin_field = y_bin_field
        self._x_log = x_log
        self._y_log = y_log
        self.end_collect = end_collect
        self.x_n_bins = x_n_bins
        self.y_n_bins = y_n_bins

        func = {True:na.logspace, False:na.linspace}[x_log]
        bounds = fix_bounds(x_lower_bound, x_upper_bound, x_log)
        self._x_bins = func(bounds[0], bounds[1], x_n_bins + 1)
        self[x_bin_field] = self._x_bins

        func = {True:na.logspace, False:na.linspace}[y_log]
        bounds = fix_bounds(y_lower_bound, y_upper_bound, y_log)
        self._y_bins = func(bounds[0], bounds[1], y_n_bins + 1)
        self[y_bin_field] = self._y_bins

        if na.any(na.isnan(self[x_bin_field])) \
            or na.any(na.isnan(self[y_bin_field])):
            mylog.error("Your min/max values for x, y have given me a nan.")
            mylog.error("Usually this means you are asking for log, with a zero bound.")
            raise ValueError

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
                binned_field = na.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = na.add.accumulate(binned_field, axis=1)
        return binned_field, weight_field, q_field, \
            used_field.astype("bool")

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data_x = self._get_field(source, self.x_bin_field, check_cut)
        source_data_y = self._get_field(source, self.y_bin_field, check_cut)
        if source_data_x.size == 0:
            raise EmptyProfileData()

        if self.end_collect:
            mi = na.arange(source_data_x.size)
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
        if self.end_collect:
            bin_indices_x = na.minimum(na.maximum(1, bin_indices_x), self.x_n_bins) - 1
            bin_indices_y = na.minimum(na.maximum(1, bin_indices_y), self.y_n_bins) - 1

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
            if self._x_log: x=na.log10(x)
            if self._y_log: y=na.log10(y)
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
        x,y = na.meshgrid(x,y)
        field_data = [x.ravel(), y.ravel()]
        if bin_style is not 'both':
            field_data += [self.field_data[field][:-1,:-1].ravel() for field in fields
                           if field not in [self.x_bin_field, self.y_bin_field]]
        else:
            field_data += [self.field_data[field].ravel() for field in fields
                           if field not in [self.x_bin_field, self.y_bin_field]]

        field_data = na.array(field_data)
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
    if logit: return na.log10(upper), na.log10(lower)
    return upper, lower

class BinnedProfile2DInlineCut(BinnedProfile2D):
    def __init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 end_collect=False):
        self.indices = data_source["Ones"].astype("bool")
        BinnedProfile2D.__init__(self, data_source,
                 x_n_bins, x_bin_field, x_lower_bound, x_upper_bound, x_log,
                 y_n_bins, y_bin_field, y_lower_bound, y_upper_bound, y_log,
                 end_collect)

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
        bin_indices_x = args[1][self.indices].ravel().astype('int64')
        bin_indices_y = args[2][self.indices].ravel().astype('int64')
        source_data = source_data[mi][self.indices]
        weight_data = weight_data[mi][self.indices]
        nx = bin_indices_x.size
        #mylog.debug("Binning %s / %s times", source_data.size, nx)
        Bin2DProfile(bin_indices_x, bin_indices_y, weight_data, source_data,
                     weight_field, binned_field, used_field)
        if accumulation: # Fix for laziness
            if not iterable(accumulation):
                raise SyntaxError("Accumulation needs to have length 2")
            if accumulation[0]:
                binned_field = na.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = na.add.accumulate(binned_field, axis=1)
        return binned_field, weight_field, used_field.astype('bool')

        
class BinnedProfile3D(BinnedProfile):
    """
    A 'Profile' produces either a weighted (or unweighted) average
    or a straight sum of a field in a bin defined by two other
    fields.  In the case of a weighted average, we have: p_i =
    sum( w_i * v_i ) / sum(w_i)
    
    We accept a *data_source*, which will be binned into
    *(x,y,z)_n_bins* by the field *(x,y,z)_bin_field* between the
    *(x,y,z)_lower_bound* and the *(x,y,z)_upper_bound*.  These bins may or
    may not be equally divided in log-space as specified by
    *(x,y,z)_log*, and the *lazy_reader* flag controls
    whether we use a memory conservative approach. If
    *end_collect* is True, take all values outside the given
    bounds and store them in the 0 and *n_bins*-1 values.
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

        func = {True:na.logspace, False:na.linspace}[x_log]
        bounds = fix_bounds(x_lower_bound, x_upper_bound, x_log)
        self._x_bins = func(bounds[0], bounds[1], x_n_bins + 1)
        self[x_bin_field] = self._x_bins

        func = {True:na.logspace, False:na.linspace}[y_log]
        bounds = fix_bounds(y_lower_bound, y_upper_bound, y_log)
        self._y_bins = func(bounds[0], bounds[1], y_n_bins + 1)
        self[y_bin_field] = self._y_bins

        func = {True:na.logspace, False:na.linspace}[z_log]
        bounds = fix_bounds(z_lower_bound, z_upper_bound, z_log)
        self._z_bins = func(bounds[0], bounds[1], z_n_bins + 1)
        self[z_bin_field] = self._z_bins

        if na.any(na.isnan(self[x_bin_field])) \
            or na.any(na.isnan(self[y_bin_field])) \
            or na.any(na.isnan(self[z_bin_field])):
            mylog.error("Your min/max values for x, y or z have given me a nan.")
            mylog.error("Usually this means you are asking for log, with a zero bound.")
            raise ValueError

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
                binned_field = na.add.accumulate(binned_field, axis=0)
            if accumulation[1]:
                binned_field = na.add.accumulate(binned_field, axis=1)
            if accumulation[2]:
                binned_field = na.add.accumulate(binned_field, axis=2)
        return binned_field, weight_field, q_field, \
            used_field.astype("bool")

    @preserve_source_parameters
    def _get_bins(self, source, check_cut=False):
        source_data_x = self._get_field(source, self.x_bin_field, check_cut)
        source_data_y = self._get_field(source, self.y_bin_field, check_cut)
        source_data_z = self._get_field(source, self.z_bin_field, check_cut)
        if source_data_x.size == 0:
            raise EmptyProfileData()
        if self.end_collect:
            mi = na.arange(source_data_x.size)
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
            raise EmptyProfileData()

        bin_indices_x = na.digitize(sd_x, self._x_bins) - 1
        bin_indices_y = na.digitize(sd_y, self._y_bins) - 1
        bin_indices_z = na.digitize(sd_z, self._z_bins) - 1
        if self.end_collect:
            bin_indices_x = na.minimum(na.maximum(1, bin_indices_x), self.x_n_bins) - 1
            bin_indices_y = na.minimum(na.maximum(1, bin_indices_y), self.y_n_bins) - 1
            bin_indices_z = na.minimum(na.maximum(1, bin_indices_z), self.z_n_bins) - 1

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
            if self._x_log: x=na.log10(x)
            if self._y_log: y=na.log10(y)
            if self._z_log: z=na.log10(z)
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
        for field in self.field_data:
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
        self.field_data = YTFieldData()
        prof_arr = pf.h.get_data("/Profiles", name)
        if prof_arr is None: raise KeyError("No such array")
        for ax in 'xyz':
            for base in ['%s_bin_field', '_%s_log']:
                setattr(self, base % ax, prof_arr.getAttr(base % ax))
        for ax in 'xyz':
            fn = getattr(self, '%s_bin_field' % ax)
            self.field_data[fn] = prof_arr.getAttr('%s_bin_values' % ax)
        shape = prof_arr.getAttr('shape')
        for fn, fd in zip(prof_arr.getAttr('field_order'),
                          prof_arr.read().transpose()):
            self.field_data[fn] = fd.reshape(shape)

    def add_fields(self, *args, **kwargs):
        raise RuntimeError("Sorry, you can't add to a stored profile.")

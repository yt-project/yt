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

import numpy as np

from yt.fields.derived_field import DerivedField
from yt.frontends.ytdata.utilities import \
    save_as_dataset
from yt.funcs import \
    get_output_filename, \
    ensure_list, \
    iterable
from yt.units.yt_array import \
    array_like_field, \
    YTQuantity
from yt.units.unit_object import Unit
from yt.data_objects.data_containers import YTFieldData
from yt.utilities.exceptions import \
    YTIllDefinedProfile
from yt.utilities.lib.misc_utilities import \
    new_bin_profile1d, \
    new_bin_profile2d, \
    new_bin_profile3d
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_objects
from yt.utilities.lib.particle_mesh_operations import \
    CICDeposit_2, \
    NGPDeposit_2


def _sanitize_min_max_units(amin, amax, finfo, registry):
    # returns a copy of amin and amax, converted to finfo's output units
    umin = getattr(amin, 'units', None)
    umax = getattr(amax, 'units', None)
    if umin is None:
        umin = Unit(finfo.output_units, registry=registry)
        rmin = YTQuantity(amin, umin)
    else:
        rmin = amin.in_units(finfo.output_units)
    if umax is None:
        umax = Unit(finfo.output_units, registry=registry)
        rmax = YTQuantity(amax, umax)
    else:
        rmax = amax.in_units(finfo.output_units)
    return rmin, rmax

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
        self.field_info = {}
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
        for f in fields:
            self.field_info[f] = self.data_source.ds.field_info[f]
        temp_storage = ProfileFieldAccumulator(len(fields), self.size)
        citer = self.data_source.chunks([], "io")
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
            units = chunk.ds.field_info[field].output_units
            arr[:,i] = chunk[field][filter].in_units(units)
        if self.weight_field is not None:
            units = chunk.ds.field_info[self.weight_field].output_units
            weight_data = chunk[self.weight_field].in_units(units)
        else:
            weight_data = np.ones(filter.size, dtype="float64")
        weight_data = weight_data[filter]
        # So that we can pass these into
        return arr, weight_data, bin_fields

    def __getitem__(self, field):
        fname = self.field_map.get(field, None)
        if fname is None:
            if isinstance(field, tuple):
                fname = self.field_map.get(field[1], None)
            elif isinstance(field, DerivedField):
                fname = self.field_map.get(field.name[1], None)
        if fname is None:
            raise KeyError(field)
        else:
            if getattr(self, 'fractional', False):
                return self.field_data[fname]
            else:
                return self.field_data[fname].in_units(self.field_units[fname])

    def items(self):
        return [(k,self[k]) for k in self.field_data.keys()]

    def keys(self):
        return self.field_data.keys()

    def __iter__(self):
        return sorted(self.items())

    def _get_bins(self, mi, ma, n, take_log):
        if take_log:
            ret = np.logspace(np.log10(mi), np.log10(ma), n+1)
            # at this point ret[0] and ret[-1] are not exactly equal to
            # mi and ma due to round-off error. Let's force them to be
            # mi and ma exactly to avoid incorrectly discarding cells near
            # the edges. See Issue #1300.
            ret[0], ret[-1] = mi, ma
            return ret
        else:
            return np.linspace(mi, ma, n+1)

    def save_as_dataset(self, filename=None):
        r"""Export a profile to a reloadable yt dataset.

        This function will take a profile and output a dataset
        containing all relevant fields.  The resulting dataset
        can be reloaded as a yt dataset.

        Parameters
        ----------
        filename : str, optional
            The name of the file to be written.  If None, the name
            will be a combination of the original dataset plus
            the type of object, e.g., Profile1D.

        Returns
        -------
        filename : str
            The name of the file that has been created.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
        >>> ad = ds.all_data()
        >>> profile = yt.create_profile(ad, ["density", "temperature"],
        ...                            "cell_mass", weight_field=None,
        ...                             n_bins=(128, 128))
        >>> fn = profile.save_as_dataset()
        >>> prof_ds = yt.load(fn)
        >>> print (prof_ds.data["cell_mass"])
        (128, 128)
        >>> print (prof_ds.data["x"].shape) # x bins as 1D array
        (128,)
        >>> print (prof_ds.data["density"]) # x bins as 2D array
        (128, 128)
        >>> p = yt.PhasePlot(prof_ds.data, "density", "temperature",
        ...                  "cell_mass", weight_field=None)
        >>> p.save()

        """

        keyword = "%s_%s" % (str(self.ds), self.__class__.__name__)
        filename = get_output_filename(filename, keyword, ".h5")

        args = ("field", "log")
        extra_attrs = {"data_type": "yt_profile",
                       "profile_dimensions": self.size,
                       "weight_field": self.weight_field,
                       "fractional": self.fractional}
        data = {}
        data.update(self.field_data)
        data["weight"] = self.weight
        data["used"] = self.used.astype("float64")

        dimensionality = 0
        bin_data = []
        for ax in "xyz":
            if hasattr(self, ax):
                dimensionality += 1
                data[ax] = getattr(self, ax)
                bin_data.append(data[ax])
                bin_field_name = "%s_bins" % ax
                data[bin_field_name] = getattr(self, bin_field_name)
                extra_attrs["%s_range" % ax] = self.ds.arr([data[bin_field_name][0],
                                                            data[bin_field_name][-1]])
                for arg in args:
                    key = "%s_%s" % (ax, arg)
                    extra_attrs[key] = getattr(self, key)

        bin_fields = np.meshgrid(*bin_data)
        for i, ax in enumerate("xyz"[:dimensionality]):
            data[getattr(self, "%s_field" % ax)] = bin_fields[i]

        extra_attrs["dimensionality"] = dimensionality
        ftypes = dict([(field, "data") for field in data])
        save_as_dataset(self.ds, filename, data, field_types=ftypes,
                        extra_attrs=extra_attrs)

        return filename

class ProfileNDFromDataset(ProfileND):
    """
    An ND profile object loaded from a ytdata dataset.
    """
    def __init__(self, ds):
        ProfileND.__init__(self, ds.data, ds.parameters["weight_field"])
        self.fractional = ds.parameters["fractional"]
        exclude_fields = ["used", "weight"]
        for ax in "xyz"[:ds.dimensionality]:
            setattr(self, ax, ds.data[ax])
            setattr(self, "%s_bins" % ax, ds.data["%s_bins" % ax])
            field_name = tuple(ds.parameters["%s_field" % ax])
            setattr(self, "%s_field" % ax, field_name)
            self.field_info[field_name] = ds.field_info[field_name]
            setattr(self, "%s_log" % ax, ds.parameters["%s_log" % ax])
            exclude_fields.extend([ax, "%s_bins" % ax,
                                   ds.parameters["%s_field" % ax][1]])
        self.weight = ds.data["weight"]
        self.used = ds.data["used"].d.astype(bool)
        profile_fields = [f for f in ds.field_list
                          if f[1] not in exclude_fields]
        for field in profile_fields:
            self.field_map[field[1]] = field
            self.field_data[field] = ds.data[field]
            self.field_info[field] = ds.field_info[field]
            self.field_units[field] = ds.data[field].units

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
        The minimum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_max : float
        The maximum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_log : boolean
        Controls whether or not the bins for the x field are evenly
        spaced in linear (False) or log (True) space.
    weight_field : string field name
        The field to weight the profiled fields by.

    """
    def __init__(self, data_source, x_field, x_n, x_min, x_max, x_log,
                 weight_field = None):
        super(Profile1D, self).__init__(data_source, weight_field)
        self.x_field = data_source._determine_fields(x_field)[0]
        self.field_info[self.x_field] = \
            self.data_source.ds.field_info[self.x_field]
        self.x_log = x_log
        x_min, x_max = _sanitize_min_max_units(
            x_min, x_max, self.field_info[self.x_field], self.ds.unit_registry)
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
        bf_x.convert_to_units(self.field_info[self.x_field].output_units)
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

class Profile1DFromDataset(ProfileNDFromDataset, Profile1D):
    """
    A 1D profile object loaded from a ytdata dataset.
    """

    def __init(self, ds):
        ProfileNDFromDataset.__init__(self, ds)

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
        The minimum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_max : float
        The maximum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_log : boolean
        Controls whether or not the bins for the x field are evenly
        spaced in linear (False) or log (True) space.
    y_field : string field name
        The field to profile as a function of along the y axis
    y_n : integer
        The number of bins along the y direction.
    y_min : float
        The minimum value of the y profile field. If supplied without units,
        assumed to be in the output units for y_field.
    y_max : float
        The maximum value of the y profile field. If supplied without units,
        assumed to be in the output units for y_field.
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
        self.x_field = data_source._determine_fields(x_field)[0]
        self.x_log = x_log
        self.field_info[self.x_field] = \
            self.data_source.ds.field_info[self.x_field]
        x_min, x_max = _sanitize_min_max_units(
            x_min, x_max, self.field_info[self.x_field], self.ds.unit_registry)
        self.x_bins = array_like_field(data_source,
                                       self._get_bins(x_min, x_max, x_n, x_log),
                                       self.x_field)
        # Y
        self.y_field = data_source._determine_fields(y_field)[0]
        self.y_log = y_log
        self.field_info[self.y_field] = \
            self.data_source.ds.field_info[self.y_field]
        y_min, y_max = _sanitize_min_max_units(
            y_min, y_max, self.field_info[self.y_field], self.ds.unit_registry)
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
        bf_x.convert_to_units(self.field_info[self.x_field].output_units)
        bin_ind_x = np.digitize(bf_x, self.x_bins) - 1
        bf_y.convert_to_units(self.field_info[self.y_field].output_units)
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

class Profile2DFromDataset(ProfileNDFromDataset, Profile2D):
    """
    A 2D profile object loaded from a ytdata dataset.
    """

    def __init(self, ds):
        ProfileNDFromDataset.__init__(self, ds)

class ParticleProfile(Profile2D):
    """An object that represents a *deposited* 2D profile. This is like a
    Profile2D, except that it is intended for particle data. Instead of just
    binning the particles, the added fields will be deposited onto the mesh
    using either the nearest-grid-point or cloud-in-cell interpolation kernels.

    Parameters
    ----------

    data_source : AMD3DData object
        The data object to be profiled
    x_field : string field name
        The field to profile as a function of along the x axis.
    x_n : integer
        The number of bins along the x direction.
    x_min : float
        The minimum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_max : float
        The maximum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    y_field : string field name
        The field to profile as a function of along the y axis
    y_n : integer
        The number of bins along the y direction.
    y_min : float
        The minimum value of the y profile field. If supplied without units,
        assumed to be in the output units for y_field.
    y_max : float
        The maximum value of the y profile field. If supplied without units,
        assumed to be in the output units for y_field.
    weight_field : string field name
        The field to use for weighting. Default is None.
    deposition : string, optional
        The interpolation kernal to be used for
        deposition. Valid choices:
        "ngp" : nearest grid point interpolation
        "cic" : cloud-in-cell interpolation

    """
    accumulation = False
    fractional = False

    def __init__(self, data_source,
                 x_field, x_n, x_min, x_max,
                 y_field, y_n, y_min, y_max,
                 weight_field=None, deposition="ngp"):

        x_field = data_source._determine_fields(x_field)[0]
        y_field = data_source._determine_fields(y_field)[0]

        # set the log parameters to False (since that doesn't make much sense
        # for deposited data) and also turn off the weight field.
        super(ParticleProfile, self).__init__(data_source,
                                              x_field,
                                              x_n, x_min, x_max, False,
                                              y_field,
                                              y_n, y_min, y_max, False,
                                              weight_field=weight_field)

        self.LeftEdge = [self.x_bins[0], self.y_bins[0]]
        self.dx = (self.x_bins[-1] - self.x_bins[0]) / x_n
        self.dy = (self.y_bins[-1] - self.y_bins[0]) / y_n
        self.CellSize = [self.dx, self.dy]
        self.CellVolume = np.product(self.CellSize)
        self.GridDimensions = np.array([x_n, y_n], dtype=np.int32)
        self.known_styles = ["ngp", "cic"]
        if deposition not in self.known_styles:
            raise NotImplementedError(deposition)
        self.deposition = deposition

    # Either stick the particle field in the nearest bin,
    # or spread it out using the 2D CIC deposition function
    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None: return
        fdata, wdata, (bf_x, bf_y) = rv
        # make sure everything has the same units before deposition.
        # the units will be scaled to the correct values later.
        LE = np.array([self.LeftEdge[0].in_units(bf_x.units),
                       self.LeftEdge[1].in_units(bf_y.units)])
        cell_size = np.array([self.CellSize[0].in_units(bf_x.units),
                              self.CellSize[1].in_units(bf_y.units)])
        for fi, field in enumerate(fields):
            Np = fdata[:, fi].size

            if self.deposition == "ngp":
                func = NGPDeposit_2
            elif self.deposition == 'cic':
                func = CICDeposit_2

            if self.weight_field is None:
                deposit_vals = fdata[:, fi]
            else:
                deposit_vals = wdata*fdata[:, fi]

            func(bf_x, bf_y, deposit_vals, Np,
                 storage.values[:, :, fi],
                 LE,
                 self.GridDimensions,
                 cell_size)

            locs = storage.values[:, :, fi] > 0.0
            storage.used[locs] = True

            if self.weight_field is not None:
                func(bf_x, bf_y, wdata, Np,
                     storage.weight_values,
                     LE,
                     self.GridDimensions,
                     cell_size)
            else:
                storage.weight_values[locs] = 1.0
            storage.mvalues[locs, fi] = storage.values[locs, fi] \
                                        / storage.weight_values[locs]
        # We've binned it!


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
        The minimum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_max : float
        The maximum value of the x profile field. If supplied without units,
        assumed to be in the output units for x_field.
    x_log : boolean
        Controls whether or not the bins for the x field are evenly
        spaced in linear (False) or log (True) space.
    y_field : string field name
        The field to profile as a function of along the y axis
    y_n : integer
        The number of bins along the y direction.
    y_min : float
        The minimum value of the y profile field. If supplied without units,
        assumed to be in the output units for y_field.
    y_max : float
        The maximum value of the y profile field. If supplied without units,
        assumed to be in the output units for y_field.
    y_log : boolean
        Controls whether or not the bins for the y field are evenly
        spaced in linear (False) or log (True) space.
    z_field : string field name
        The field to profile as a function of along the z axis
    z_n : integer
        The number of bins along the z direction.
    z_min : float
        The minimum value of the z profile field. If supplied without units,
        assumed to be in the output units for z_field.
    z_max : float
        The maximum value of thee z profile field. If supplied without units,
        assumed to be in the output units for z_field.
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
        self.x_field = data_source._determine_fields(x_field)[0]
        self.x_log = x_log
        self.field_info[self.x_field] = \
            self.data_source.ds.field_info[self.x_field]
        x_min, x_max = _sanitize_min_max_units(
            x_min, x_max, self.field_info[self.x_field], self.ds.unit_registry)
        self.x_bins = array_like_field(data_source,
                                       self._get_bins(x_min, x_max, x_n, x_log),
                                       self.x_field)
        # Y
        self.y_field = data_source._determine_fields(y_field)[0]
        self.y_log = y_log
        self.field_info[self.y_field] = \
            self.data_source.ds.field_info[self.y_field]
        y_min, y_max = _sanitize_min_max_units(
            y_min, y_max, self.field_info[self.y_field], self.ds.unit_registry)
        self.y_bins = array_like_field(data_source,
                                       self._get_bins(y_min, y_max, y_n, y_log),
                                       self.y_field)
        # Z
        self.z_field = data_source._determine_fields(z_field)[0]
        self.z_log = z_log
        self.field_info[self.z_field] = \
            self.data_source.ds.field_info[self.z_field]
        z_min, z_max = _sanitize_min_max_units(
            z_min, z_max, self.field_info[self.z_field], self.ds.unit_registry)
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
        bf_x.convert_to_units(self.field_info[self.x_field].output_units)
        bin_ind_x = np.digitize(bf_x, self.x_bins) - 1
        bf_y.convert_to_units(self.field_info[self.y_field].output_units)
        bin_ind_y = np.digitize(bf_y, self.y_bins) - 1
        bf_z.convert_to_units(self.field_info[self.z_field].output_units)
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

class Profile3DFromDataset(ProfileNDFromDataset, Profile3D):
    """
    A 2D profile object loaded from a ytdata dataset.
    """

    def __init(self, ds):
        ProfileNDFromDataset.__init__(self, ds)

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
                   accumulation=False, fractional=False,
                   deposition='ngp'):
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
    deposition : Controls the type of deposition used for ParticlePhasePlots.
        Valid choices are 'ngp' and 'cic'. Default is 'ngp'. This parameter is
        ignored the if the input fields are not of particle type.


    Examples
    --------

    Create a 1d profile.  Access bin field from profile.x and field
    data from profile[<field_name>].

    >>> ds = load("DD0046/DD0046")
    >>> ad = ds.h.all_data()
    >>> profile = create_profile(ad, [("gas", "density")],
    ...                              [("gas", "temperature"),
    ...                               ("gas", "velocity_x")])
    >>> print (profile.x)
    >>> print (profile["gas", "temperature"])

    """
    bin_fields = data_source._determine_fields(bin_fields)
    fields = ensure_list(fields)
    is_pfield = [data_source.ds._get_field_info(f).particle_type
                 for f in bin_fields + fields]
    wf = None
    if weight_field is not None:
        wf = data_source.ds._get_field_info(weight_field)
        is_pfield.append(wf.particle_type)
        wf = wf.name

    if any(is_pfield) and not all(is_pfield):
        raise YTIllDefinedProfile(
            bin_fields, data_source._determine_fields(fields), wf, is_pfield)
    elif len(bin_fields) == 1:
        cls = Profile1D
    elif len(bin_fields) == 2 and all(is_pfield):
        # log bin_fields set to False for Particle Profiles.
        # doesn't make much sense for CIC deposition.
        # accumulation and fractional set to False as well.
        logs = {bin_fields[0]: False, bin_fields[1]: False}
        accumulation = False
        fractional = False
        cls = ParticleProfile
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
    if weight_field is not None and cls == ParticleProfile:
        weight_field, = data_source._determine_fields([weight_field])
        if not data_source.ds._get_field_info(weight_field).particle_type:
            weight_field = None
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
        # pad extrema by epsilon so cells at bin edges are not excluded
        for i, (mi, ma) in enumerate(ex):
            mi = mi - np.spacing(mi)
            ma = ma + np.spacing(ma)
            ex[i][0], ex[i][1] = mi, ma
    else:
        ex = []
        for bin_field in bin_fields:
            bf_units = data_source.ds.field_info[bin_field].output_units
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
    if cls is ParticleProfile:
        args = [data_source]
        for f, n, (mi, ma) in zip(bin_fields, n_bins, ex):
            args += [f, n, mi, ma]
        obj = cls(*args, weight_field=weight_field, deposition=deposition)
    else:
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
        for field, unit in units.items():
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

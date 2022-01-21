import numpy as np
from more_itertools import collapse

from yt.data_objects.field_data import YTFieldData
from yt.fields.derived_field import DerivedField
from yt.frontends.ytdata.utilities import save_as_dataset
from yt.funcs import get_output_filename, is_sequence, iter_fields, mylog
from yt.units.unit_object import Unit  # type: ignore
from yt.units.yt_array import YTQuantity, array_like_field
from yt.utilities.exceptions import (
    YTIllDefinedBounds,
    YTIllDefinedProfile,
    YTProfileDataShape,
)
from yt.utilities.lib.misc_utilities import (
    new_bin_profile1d,
    new_bin_profile2d,
    new_bin_profile3d,
)
from yt.utilities.lib.particle_mesh_operations import CICDeposit_2, NGPDeposit_2
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    ParallelAnalysisInterface,
    parallel_objects,
)


def _sanitize_min_max_units(amin, amax, finfo, registry):
    # returns a copy of amin and amax, converted to finfo's output units
    umin = getattr(amin, "units", None)
    umax = getattr(amax, "units", None)
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
        if hasattr(source, "field_parameters"):
            old_params = source.field_parameters
            source.field_parameters = prof._data_source.field_parameters
            tr = func(*args, **kwargs)
            source.field_parameters = old_params
        else:
            tr = func(*args, **kwargs)
        return tr

    return save_state


class ProfileFieldAccumulator:
    def __init__(self, n_fields, size):
        shape = size + (n_fields,)
        self.values = np.zeros(shape, dtype="float64")
        self.mvalues = np.zeros(shape, dtype="float64")
        self.qvalues = np.zeros(shape, dtype="float64")
        self.used = np.zeros(size, dtype="bool")
        self.weight_values = np.zeros(size, dtype="float64")


class ProfileND(ParallelAnalysisInterface):
    """The profile object class"""

    def __init__(self, data_source, weight_field=None):
        self.data_source = data_source
        self.ds = data_source.ds
        self.field_map = {}
        self.field_info = {}
        self.field_data = YTFieldData()
        if weight_field is not None:
            self.standard_deviation = YTFieldData()
            weight_field = self.data_source._determine_fields(weight_field)[0]
        else:
            self.standard_deviation = None
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
            self.field_units[field] = Unit(new_unit, registry=self.ds.unit_registry)
        else:
            fd = self.field_map[field]
            if fd in self.field_units:
                self.field_units[fd] = Unit(new_unit, registry=self.ds.unit_registry)
            else:
                raise KeyError(f"{field} not in profile!")

    def _finalize_storage(self, fields, temp_storage):
        # We use our main comm here
        # This also will fill _field_data

        for i, _field in enumerate(fields):
            # q values are returned as q * weight but we want just q
            temp_storage.qvalues[..., i][
                temp_storage.used
            ] /= temp_storage.weight_values[temp_storage.used]

        # get the profile data from all procs
        all_store = {self.comm.rank: temp_storage}
        all_store = self.comm.par_combine_object(all_store, "join", datatype="dict")

        all_val = np.zeros_like(temp_storage.values)
        all_mean = np.zeros_like(temp_storage.mvalues)
        all_std = np.zeros_like(temp_storage.qvalues)
        all_weight = np.zeros_like(temp_storage.weight_values)
        all_used = np.zeros_like(temp_storage.used, dtype="bool")

        # Combine the weighted mean and standard deviation from each processor.
        # For two samples with total weight, mean, and standard deviation
        # given by w, m, and s, their combined mean and standard deviation are:
        # m12 = (m1 * w1 + m2 * w2) / (w1 + w2)
        # s12 = (m1 * (s1**2 + (m1 - m12)**2) +
        #        m2 * (s2**2 + (m2 - m12)**2)) / (w1 + w2)
        # Here, the mvalues are m and the qvalues are s**2.
        for p in sorted(all_store.keys()):
            all_used += all_store[p].used
            old_mean = all_mean.copy()
            old_weight = all_weight.copy()
            all_weight[all_store[p].used] += all_store[p].weight_values[
                all_store[p].used
            ]
            for i, _field in enumerate(fields):
                all_val[..., i][all_store[p].used] += all_store[p].values[..., i][
                    all_store[p].used
                ]

                all_mean[..., i][all_store[p].used] = (
                    all_mean[..., i] * old_weight
                    + all_store[p].mvalues[..., i] * all_store[p].weight_values
                )[all_store[p].used] / all_weight[all_store[p].used]

                all_std[..., i][all_store[p].used] = (
                    old_weight
                    * (all_std[..., i] + (old_mean[..., i] - all_mean[..., i]) ** 2)
                    + all_store[p].weight_values
                    * (
                        all_store[p].qvalues[..., i]
                        + (all_store[p].mvalues[..., i] - all_mean[..., i]) ** 2
                    )
                )[all_store[p].used] / all_weight[all_store[p].used]

        all_std = np.sqrt(all_std)
        del all_store
        self.used = all_used
        blank = ~all_used

        self.weight = all_weight
        self.weight[blank] = 0.0

        for i, field in enumerate(fields):
            if self.weight_field is None:
                self.field_data[field] = array_like_field(
                    self.data_source, all_val[..., i], field
                )
            else:
                self.field_data[field] = array_like_field(
                    self.data_source, all_mean[..., i], field
                )
                self.standard_deviation[field] = array_like_field(
                    self.data_source, all_std[..., i], field
                )
                self.standard_deviation[field][blank] = 0.0
                self.weight = array_like_field(
                    self.data_source, self.weight, self.weight_field
                )
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
        pfilter = np.ones(bin_fields[0].shape, dtype="bool")
        for (mi, ma), data in zip(self.bounds, bin_fields):
            pfilter &= data > mi
            pfilter &= data < ma
        return pfilter, [data[pfilter] for data in bin_fields]

    def _get_data(self, chunk, fields):
        # We are using chunks now, which will manage the field parameters and
        # the like.
        bin_fields = [chunk[bf] for bf in self.bin_fields]
        for i in range(1, len(bin_fields)):
            if bin_fields[0].shape != bin_fields[i].shape:
                raise YTProfileDataShape(
                    self.bin_fields[0],
                    bin_fields[0].shape,
                    self.bin_fields[i],
                    bin_fields[i].shape,
                )
        # We want to make sure that our fields are within the bounds of the
        # binning
        pfilter, bin_fields = self._filter(bin_fields)
        if not np.any(pfilter):
            return None
        arr = np.zeros((bin_fields[0].size, len(fields)), dtype="float64")
        for i, field in enumerate(fields):
            if pfilter.shape != chunk[field].shape:
                raise YTProfileDataShape(
                    self.bin_fields[0], bin_fields[0].shape, field, chunk[field].shape
                )
            units = chunk.ds.field_info[field].output_units
            arr[:, i] = chunk[field][pfilter].in_units(units)
        if self.weight_field is not None:
            if pfilter.shape != chunk[self.weight_field].shape:
                raise YTProfileDataShape(
                    self.bin_fields[0],
                    bin_fields[0].shape,
                    self.weight_field,
                    chunk[self.weight_field].shape,
                )
            units = chunk.ds.field_info[self.weight_field].output_units
            weight_data = chunk[self.weight_field].in_units(units)
        else:
            weight_data = np.ones(pfilter.shape, dtype="float64")
        weight_data = weight_data[pfilter]
        # So that we can pass these into
        return arr, weight_data, bin_fields

    def __getitem__(self, field):
        if field in self.field_data:
            fname = field
        else:
            # deal with string vs tuple field names and attempt to guess which field
            # we are supposed to be talking about
            fname = self.field_map.get(field, None)
            if isinstance(field, tuple):
                fname = self.field_map.get(field[1], None)
                if fname != field:
                    raise KeyError(
                        "Asked for field '{}' but only have data for "
                        "fields '{}'".format(field, list(self.field_data.keys()))
                    )
            elif isinstance(field, DerivedField):
                fname = self.field_map.get(field.name[1], None)
            if fname is None:
                raise KeyError(field)
        if getattr(self, "fractional", False):
            return self.field_data[fname]
        else:
            return self.field_data[fname].in_units(self.field_units[fname])

    def items(self):
        return [(k, self[k]) for k in self.field_data.keys()]

    def keys(self):
        return self.field_data.keys()

    def __iter__(self):
        return sorted(self.items())

    def _get_bins(self, mi, ma, n, take_log):
        if take_log:
            ret = np.logspace(np.log10(mi), np.log10(ma), n + 1)
            # at this point ret[0] and ret[-1] are not exactly equal to
            # mi and ma due to round-off error. Let's force them to be
            # mi and ma exactly to avoid incorrectly discarding cells near
            # the edges. See Issue #1300.
            ret[0], ret[-1] = mi, ma
            return ret
        else:
            return np.linspace(mi, ma, n + 1)

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
        >>> profile = yt.create_profile(
        ...     ad,
        ...     [("gas", "density"), ("gas", "temperature")],
        ...     ("gas", "mass"),
        ...     weight_field=None,
        ...     n_bins=(128, 128),
        ... )
        >>> fn = profile.save_as_dataset()
        >>> prof_ds = yt.load(fn)
        >>> print(prof_ds.data[("gas", "mass")])
        (128, 128)
        >>> print(prof_ds.data[("index", "x")].shape)  # x bins as 1D array
        (128,)
        >>> print(prof_ds.data[("gas", "density")])  # x bins as 2D array
        (128, 128)
        >>> p = yt.PhasePlot(
        ...     prof_ds.data,
        ...     ("gas", "density"),
        ...     ("gas", "temperature"),
        ...     ("gas", "mass"),
        ...     weight_field=None,
        ... )
        >>> p.save()

        """

        keyword = f"{str(self.ds)}_{self.__class__.__name__}"
        filename = get_output_filename(filename, keyword, ".h5")

        args = ("field", "log")
        extra_attrs = {
            "data_type": "yt_profile",
            "profile_dimensions": self.size,
            "weight_field": self.weight_field,
            "fractional": self.fractional,
            "accumulation": self.accumulation,
        }
        data = {}
        data.update(self.field_data)
        data["weight"] = self.weight
        data["used"] = self.used.astype("float64")
        std = "standard_deviation"
        if self.weight_field is not None:
            std_data = getattr(self, std)
            data.update({(std, field[1]): std_data[field] for field in self.field_data})

        dimensionality = 0
        bin_data = []
        for ax in "xyz":
            if hasattr(self, ax):
                dimensionality += 1
                data[ax] = getattr(self, ax)
                bin_data.append(data[ax])
                bin_field_name = f"{ax}_bins"
                data[bin_field_name] = getattr(self, bin_field_name)
                extra_attrs[f"{ax}_range"] = self.ds.arr(
                    [data[bin_field_name][0], data[bin_field_name][-1]]
                )
                for arg in args:
                    key = f"{ax}_{arg}"
                    extra_attrs[key] = getattr(self, key)

        bin_fields = np.meshgrid(*bin_data)
        for i, ax in enumerate("xyz"[:dimensionality]):
            data[getattr(self, f"{ax}_field")] = bin_fields[i]

        extra_attrs["dimensionality"] = dimensionality
        ftypes = {field: "data" for field in data if field[0] != std}
        if self.weight_field is not None:
            ftypes.update({(std, field[1]): std for field in self.field_data})
        save_as_dataset(
            self.ds, filename, data, field_types=ftypes, extra_attrs=extra_attrs
        )

        return filename


class ProfileNDFromDataset(ProfileND):
    """
    An ND profile object loaded from a ytdata dataset.
    """

    def __init__(self, ds):
        ProfileND.__init__(self, ds.data, ds.parameters.get("weight_field", None))
        self.fractional = ds.parameters.get("fractional", False)
        self.accumulation = ds.parameters.get("accumulation", False)
        exclude_fields = ["used", "weight"]
        for ax in "xyz"[: ds.dimensionality]:
            setattr(self, ax, ds.data[("data", ax)])
            ax_bins = f"{ax}_bins"
            ax_field = f"{ax}_field"
            ax_log = f"{ax}_log"
            setattr(self, ax_bins, ds.data[("data", ax_bins)])
            field_name = tuple(ds.parameters.get(ax_field, (None, None)))
            setattr(self, ax_field, field_name)
            self.field_info[field_name] = ds.field_info[field_name]
            setattr(self, ax_log, ds.parameters.get(ax_log, False))
            exclude_fields.extend([ax, ax_bins, field_name[1]])
        self.weight = ds.data[("data", "weight")]
        self.used = ds.data[("data", "used")].d.astype(bool)
        profile_fields = [
            f
            for f in ds.field_list
            if f[1] not in exclude_fields and f[0] != "standard_deviation"
        ]
        for field in profile_fields:
            self.field_map[field[1]] = field
            self.field_data[field] = ds.data[field]
            self.field_info[field] = ds.field_info[field]
            self.field_units[field] = ds.data[field].units
            if ("standard_deviation", field[1]) in ds.field_list:
                self.standard_deviation[field] = ds.data["standard_deviation", field[1]]


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
    override_bins_x : array
        Array to set as xbins and ignore other parameters if set

    """

    def __init__(
        self,
        data_source,
        x_field,
        x_n,
        x_min,
        x_max,
        x_log,
        weight_field=None,
        override_bins_x=None,
    ):
        super().__init__(data_source, weight_field)
        self.x_field = data_source._determine_fields(x_field)[0]
        self.field_info[self.x_field] = self.data_source.ds.field_info[self.x_field]
        self.x_log = x_log
        x_min, x_max = _sanitize_min_max_units(
            x_min, x_max, self.field_info[self.x_field], self.ds.unit_registry
        )
        self.x_bins = array_like_field(
            data_source, self._get_bins(x_min, x_max, x_n, x_log), self.x_field
        )

        if override_bins_x is not None:
            self.x_bins = array_like_field(data_source, override_bins_x, self.x_field)

        self.size = (self.x_bins.size - 1,)
        self.bin_fields = (self.x_field,)
        self.x = 0.5 * (self.x_bins[1:] + self.x_bins[:-1])

    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None:
            return
        fdata, wdata, (bf_x,) = rv
        bf_x.convert_to_units(self.field_info[self.x_field].output_units)
        bin_ind = np.digitize(bf_x, self.x_bins) - 1
        new_bin_profile1d(
            bin_ind,
            wdata,
            fdata,
            storage.weight_values,
            storage.values,
            storage.mvalues,
            storage.qvalues,
            storage.used,
        )

        # We've binned it!

    def set_x_unit(self, new_unit):
        """Sets a new unit for the x field

        Parameters
        ----------
        new_unit : string or Unit object
            The name of the new unit.
        """
        self.x_bins.convert_to_units(new_unit)
        self.x = 0.5 * (self.x_bins[1:] + self.x_bins[:-1])

    @property
    def bounds(self):
        return ((self.x_bins[0], self.x_bins[-1]),)

    def plot(self):
        r"""
        This returns a :class:`~yt.visualization.profile_plotter.ProfilePlot`
        with the fields that have been added to this object.
        """
        from yt.visualization.profile_plotter import ProfilePlot

        return ProfilePlot.from_profiles(self)

    def _export_prep(self, fields, only_used):
        if only_used:
            idxs = self.used
        else:
            idxs = slice(None, None, None)
        if not only_used and not np.all(self.used):
            masked = True
        else:
            masked = False
        if fields is None:
            fields = self.field_data.keys()
        else:
            fields = self.data_source._determine_fields(fields)
        return idxs, masked, fields

    def to_dataframe(self, fields=None, only_used=False):
        r"""Export a profile object to a pandas DataFrame.

        This function will take a data object and construct from it and
        optionally a list of fields a pandas DataFrame object.  If pandas is
        not importable, this will raise ImportError.

        Parameters
        ----------
        fields : list of strings or tuple field names, default None
            If this is supplied, it is the list of fields to be exported into
            the DataFrame. If not supplied, whatever fields exist in the
            profile, along with the bin field, will be exported.
        only_used : boolean, default False
            If True, only the bins which have data will be exported. If False,
            all of the bins will be exported, but the elements for those bins
            in the data arrays will be filled with NaNs.

        Returns
        -------
        df : :class:`~pandas.DataFrame`
            The data contained in the profile.

        Examples
        --------
        >>> sp = ds.sphere("c", (0.1, "unitary"))
        >>> p = sp.profile(
        ...     ("index", "radius"), [("gas", "density"), ("gas", "temperature")]
        ... )
        >>> df1 = p.to_dataframe()
        >>> df2 = p.to_dataframe(fields=("gas", "density"), only_used=True)
        """
        from yt.utilities.on_demand_imports import _pandas as pd

        idxs, masked, fields = self._export_prep(fields, only_used)
        pdata = {self.x_field[-1]: self.x[idxs]}
        for field in fields:
            pdata[field[-1]] = self[field][idxs]
        df = pd.DataFrame(pdata)
        if masked:
            mask = np.zeros(df.shape, dtype="bool")
            mask[~self.used, 1:] = True
            df.mask(mask, inplace=True)
        return df

    def to_astropy_table(self, fields=None, only_used=False):
        """
        Export the profile data to a :class:`~astropy.table.table.QTable`,
        which is a Table object which is unit-aware. The QTable can then
        be exported to an ASCII file, FITS file, etc.

        See the AstroPy Table docs for more details:
        http://docs.astropy.org/en/stable/table/

        Parameters
        ----------
        fields : list of strings or tuple field names, default None
            If this is supplied, it is the list of fields to be exported into
            the DataFrame. If not supplied, whatever fields exist in the
            profile, along with the bin field, will be exported.
        only_used : boolean, optional
            If True, only the bins which are used are copied
            to the QTable as rows. If False, all bins are
            copied, but the bins which are not used are masked.
            Default: False

        Returns
        -------
        df : :class:`~astropy.table.QTable`
            The data contained in the profile.

        Examples
        --------
        >>> sp = ds.sphere("c", (0.1, "unitary"))
        >>> p = sp.profile(
        ...     ("index", "radius"), [("gas", "density"), ("gas", "temperature")]
        ... )
        >>> qt1 = p.to_astropy_table()
        >>> qt2 = p.to_astropy_table(fields=("gas", "density"), only_used=True)
        """
        from astropy.table import QTable

        idxs, masked, fields = self._export_prep(fields, only_used)
        qt = QTable(masked=masked)
        qt[self.x_field[-1]] = self.x[idxs].to_astropy()
        if masked:
            qt[self.x_field[-1]].mask = self.used
        for field in fields:
            qt[field[-1]] = self[field][idxs].to_astropy()
            if masked:
                qt[field[-1]].mask = self.used
        return qt


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
    override_bins_x : array
        Array to set as xbins and ignore other parameters if set
    override_bins_y : array
        Array to set as ybins and ignore other parameters if set

    """

    def __init__(
        self,
        data_source,
        x_field,
        x_n,
        x_min,
        x_max,
        x_log,
        y_field,
        y_n,
        y_min,
        y_max,
        y_log,
        weight_field=None,
        override_bins_x=None,
        override_bins_y=None,
    ):
        super().__init__(data_source, weight_field)
        # X
        self.x_field = data_source._determine_fields(x_field)[0]
        self.x_log = x_log
        self.field_info[self.x_field] = self.data_source.ds.field_info[self.x_field]
        x_min, x_max = _sanitize_min_max_units(
            x_min, x_max, self.field_info[self.x_field], self.ds.unit_registry
        )
        self.x_bins = array_like_field(
            data_source, self._get_bins(x_min, x_max, x_n, x_log), self.x_field
        )
        if override_bins_x is not None:
            self.x_bins = array_like_field(data_source, override_bins_x, self.x_field)

        # Y
        self.y_field = data_source._determine_fields(y_field)[0]
        self.y_log = y_log
        self.field_info[self.y_field] = self.data_source.ds.field_info[self.y_field]
        y_min, y_max = _sanitize_min_max_units(
            y_min, y_max, self.field_info[self.y_field], self.ds.unit_registry
        )
        self.y_bins = array_like_field(
            data_source, self._get_bins(y_min, y_max, y_n, y_log), self.y_field
        )
        if override_bins_y is not None:
            self.y_bins = array_like_field(data_source, override_bins_y, self.y_field)

        self.size = (self.x_bins.size - 1, self.y_bins.size - 1)

        self.bin_fields = (self.x_field, self.y_field)
        self.x = 0.5 * (self.x_bins[1:] + self.x_bins[:-1])
        self.y = 0.5 * (self.y_bins[1:] + self.y_bins[:-1])

    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None:
            return
        fdata, wdata, (bf_x, bf_y) = rv
        bf_x.convert_to_units(self.field_info[self.x_field].output_units)
        bin_ind_x = np.digitize(bf_x, self.x_bins) - 1
        bf_y.convert_to_units(self.field_info[self.y_field].output_units)
        bin_ind_y = np.digitize(bf_y, self.y_bins) - 1
        new_bin_profile2d(
            bin_ind_x,
            bin_ind_y,
            wdata,
            fdata,
            storage.weight_values,
            storage.values,
            storage.mvalues,
            storage.qvalues,
            storage.used,
        )
        # We've binned it!

    def set_x_unit(self, new_unit):
        """Sets a new unit for the x field

        Parameters
        ----------
        new_unit : string or Unit object
            The name of the new unit.
        """
        self.x_bins.convert_to_units(new_unit)
        self.x = 0.5 * (self.x_bins[1:] + self.x_bins[:-1])

    def set_y_unit(self, new_unit):
        """Sets a new unit for the y field

        Parameters
        ----------
        new_unit : string or Unit object
            The name of the new unit.
        """
        self.y_bins.convert_to_units(new_unit)
        self.y = 0.5 * (self.y_bins[1:] + self.y_bins[:-1])

    @property
    def bounds(self):
        return ((self.x_bins[0], self.x_bins[-1]), (self.y_bins[0], self.y_bins[-1]))

    def plot(self):
        r"""
        This returns a :class:~yt.visualization.profile_plotter.PhasePlot with
        the fields that have been added to this object.
        """
        from yt.visualization.profile_plotter import PhasePlot

        return PhasePlot.from_profile(self)


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
        The interpolation kernel to be used for
        deposition. Valid choices:
        "ngp" : nearest grid point interpolation
        "cic" : cloud-in-cell interpolation

    """

    accumulation = False
    fractional = False

    def __init__(
        self,
        data_source,
        x_field,
        x_n,
        x_min,
        x_max,
        x_log,
        y_field,
        y_n,
        y_min,
        y_max,
        y_log,
        weight_field=None,
        deposition="ngp",
    ):

        x_field = data_source._determine_fields(x_field)[0]
        y_field = data_source._determine_fields(y_field)[0]

        if deposition not in ["ngp", "cic"]:
            raise NotImplementedError(deposition)
        elif (x_log or y_log) and deposition != "ngp":
            mylog.warning(
                "cic deposition is only supported for linear axis "
                "scales, falling back to ngp deposition"
            )
            deposition = "ngp"

        self.deposition = deposition

        # set the log parameters to False (since that doesn't make much sense
        # for deposited data) and also turn off the weight field.
        super().__init__(
            data_source,
            x_field,
            x_n,
            x_min,
            x_max,
            x_log,
            y_field,
            y_n,
            y_min,
            y_max,
            y_log,
            weight_field=weight_field,
        )

    # Either stick the particle field in the nearest bin,
    # or spread it out using the 2D CIC deposition function
    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None:
            return
        fdata, wdata, (bf_x, bf_y) = rv
        # make sure everything has the same units before deposition.
        # the units will be scaled to the correct values later.

        if self.deposition == "ngp":
            func = NGPDeposit_2
        elif self.deposition == "cic":
            func = CICDeposit_2

        for fi, _field in enumerate(fields):
            if self.weight_field is None:
                deposit_vals = fdata[:, fi]
            else:
                deposit_vals = wdata * fdata[:, fi]

            field_mask = np.zeros(self.size, dtype="uint8")

            func(
                bf_x,
                bf_y,
                deposit_vals,
                fdata[:, fi].size,
                storage.values[:, :, fi],
                field_mask,
                self.x_bins,
                self.y_bins,
            )

            locs = field_mask > 0
            storage.used[locs] = True

            if self.weight_field is not None:
                func(
                    bf_x,
                    bf_y,
                    wdata,
                    fdata[:, fi].size,
                    storage.weight_values,
                    field_mask,
                    self.x_bins,
                    self.y_bins,
                )
            else:
                storage.weight_values[locs] = 1.0
            storage.mvalues[locs, fi] = (
                storage.values[locs, fi] / storage.weight_values[locs]
            )
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
    override_bins_x : array
        Array to set as xbins and ignore other parameters if set
    override_bins_y : array
        Array to set as xbins and ignore other parameters if set
    override_bins_z : array
        Array to set as xbins and ignore other parameters if set

    """

    def __init__(
        self,
        data_source,
        x_field,
        x_n,
        x_min,
        x_max,
        x_log,
        y_field,
        y_n,
        y_min,
        y_max,
        y_log,
        z_field,
        z_n,
        z_min,
        z_max,
        z_log,
        weight_field=None,
        override_bins_x=None,
        override_bins_y=None,
        override_bins_z=None,
    ):
        super().__init__(data_source, weight_field)
        # X
        self.x_field = data_source._determine_fields(x_field)[0]
        self.x_log = x_log
        self.field_info[self.x_field] = self.data_source.ds.field_info[self.x_field]
        x_min, x_max = _sanitize_min_max_units(
            x_min, x_max, self.field_info[self.x_field], self.ds.unit_registry
        )
        self.x_bins = array_like_field(
            data_source, self._get_bins(x_min, x_max, x_n, x_log), self.x_field
        )
        if override_bins_x is not None:
            self.x_bins = array_like_field(data_source, override_bins_x, self.x_field)
        # Y
        self.y_field = data_source._determine_fields(y_field)[0]
        self.y_log = y_log
        self.field_info[self.y_field] = self.data_source.ds.field_info[self.y_field]
        y_min, y_max = _sanitize_min_max_units(
            y_min, y_max, self.field_info[self.y_field], self.ds.unit_registry
        )
        self.y_bins = array_like_field(
            data_source, self._get_bins(y_min, y_max, y_n, y_log), self.y_field
        )
        if override_bins_y is not None:
            self.y_bins = array_like_field(data_source, override_bins_y, self.y_field)
        # Z
        self.z_field = data_source._determine_fields(z_field)[0]
        self.z_log = z_log
        self.field_info[self.z_field] = self.data_source.ds.field_info[self.z_field]
        z_min, z_max = _sanitize_min_max_units(
            z_min, z_max, self.field_info[self.z_field], self.ds.unit_registry
        )
        self.z_bins = array_like_field(
            data_source, self._get_bins(z_min, z_max, z_n, z_log), self.z_field
        )
        if override_bins_z is not None:
            self.z_bins = array_like_field(data_source, override_bins_z, self.z_field)

        self.size = (self.x_bins.size - 1, self.y_bins.size - 1, self.z_bins.size - 1)

        self.bin_fields = (self.x_field, self.y_field, self.z_field)
        self.x = 0.5 * (self.x_bins[1:] + self.x_bins[:-1])
        self.y = 0.5 * (self.y_bins[1:] + self.y_bins[:-1])
        self.z = 0.5 * (self.z_bins[1:] + self.z_bins[:-1])

    def _bin_chunk(self, chunk, fields, storage):
        rv = self._get_data(chunk, fields)
        if rv is None:
            return
        fdata, wdata, (bf_x, bf_y, bf_z) = rv
        bf_x.convert_to_units(self.field_info[self.x_field].output_units)
        bin_ind_x = np.digitize(bf_x, self.x_bins) - 1
        bf_y.convert_to_units(self.field_info[self.y_field].output_units)
        bin_ind_y = np.digitize(bf_y, self.y_bins) - 1
        bf_z.convert_to_units(self.field_info[self.z_field].output_units)
        bin_ind_z = np.digitize(bf_z, self.z_bins) - 1
        new_bin_profile3d(
            bin_ind_x,
            bin_ind_y,
            bin_ind_z,
            wdata,
            fdata,
            storage.weight_values,
            storage.values,
            storage.mvalues,
            storage.qvalues,
            storage.used,
        )
        # We've binned it!

    @property
    def bounds(self):
        return (
            (self.x_bins[0], self.x_bins[-1]),
            (self.y_bins[0], self.y_bins[-1]),
            (self.z_bins[0], self.z_bins[-1]),
        )

    def set_x_unit(self, new_unit):
        """Sets a new unit for the x field

        Parameters
        ----------
        new_unit : string or Unit object
            The name of the new unit.
        """
        self.x_bins.convert_to_units(new_unit)
        self.x = 0.5 * (self.x_bins[1:] + self.x_bins[:-1])

    def set_y_unit(self, new_unit):
        """Sets a new unit for the y field

        Parameters
        ----------
        new_unit : string or Unit object
            The name of the new unit.
        """
        self.y_bins.convert_to_units(new_unit)
        self.y = 0.5 * (self.y_bins[1:] + self.y_bins[:-1])

    def set_z_unit(self, new_unit):
        """Sets a new unit for the z field

        Parameters
        ----------
        new_unit : string or Unit object
            The name of the new unit.
        """
        self.z_bins.convert_to_units(new_unit)
        self.z = 0.5 * (self.z_bins[1:] + self.z_bins[:-1])


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


def create_profile(
    data_source,
    bin_fields,
    fields,
    n_bins=64,
    extrema=None,
    logs=None,
    units=None,
    weight_field=("gas", "mass"),
    accumulation=False,
    fractional=False,
    deposition="ngp",
    override_bins=None,
):
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
        each bin. Defaults to ("gas", "mass").
    accumulation : bool or list of bools
        If True, the profile values for a bin n are the cumulative sum of
        all the values from bin 0 to n.  If -True, the sum is reversed so
        that the value for bin n is the cumulative sum from bin N (total bins)
        to n.  If the profile is 2D or 3D, a list of values can be given to
        control the summation in each dimension independently.
        Default: False.
    fractional : bool
        If True the profile values are divided by the sum of all
        the profile data such that the profile represents a probability
        distribution function.
    deposition : strings
        Controls the type of deposition used for ParticlePhasePlots.
        Valid choices are 'ngp' and 'cic'. Default is 'ngp'. This parameter is
        ignored the if the input fields are not of particle type.
    override_bins : dict of bins to profile plot with
        If set, ignores n_bins and extrema settings and uses the
        supplied bins to profile the field. If a units dict is provided,
        bins are understood to be in the units specified in the dictionary.


    Examples
    --------

    Create a 1d profile.  Access bin field from profile.x and field
    data from profile[<field_name>].

    >>> ds = load("DD0046/DD0046")
    >>> ad = ds.all_data()
    >>> profile = create_profile(
    ...     ad, [("gas", "density")], [("gas", "temperature"), ("gas", "velocity_x")]
    ... )
    >>> print(profile.x)
    >>> print(profile["gas", "temperature"])

    """
    bin_fields = data_source._determine_fields(bin_fields)
    fields = list(iter_fields(fields))
    is_pfield = [
        data_source.ds._get_field_info(f).sampling_type == "particle"
        for f in bin_fields + fields
    ]
    wf = None
    if weight_field is not None:
        wf = data_source.ds._get_field_info(weight_field)
        is_pfield.append(wf.sampling_type == "particle")
        wf = wf.name

    if len(bin_fields) > 1 and isinstance(accumulation, bool):
        accumulation = [accumulation for _ in range(len(bin_fields))]

    bin_fields = data_source._determine_fields(bin_fields)
    fields = data_source._determine_fields(fields)
    units = sanitize_field_tuple_keys(units, data_source)
    extrema = sanitize_field_tuple_keys(extrema, data_source)
    logs = sanitize_field_tuple_keys(logs, data_source)
    override_bins = sanitize_field_tuple_keys(override_bins, data_source)

    if any(is_pfield) and not all(is_pfield):
        if hasattr(data_source.ds, "_sph_ptypes"):
            is_local = [
                data_source.ds.field_info[f].sampling_type == "local"
                for f in bin_fields + fields
            ]
            is_local_or_pfield = [pf or lf for (pf, lf) in zip(is_pfield, is_local)]
            if not all(is_local_or_pfield):
                raise YTIllDefinedProfile(
                    bin_fields, data_source._determine_fields(fields), wf, is_pfield
                )
        else:
            raise YTIllDefinedProfile(
                bin_fields, data_source._determine_fields(fields), wf, is_pfield
            )
    if len(bin_fields) == 1:
        cls = Profile1D
    elif len(bin_fields) == 2 and all(is_pfield):
        if deposition == "cic":
            if logs is not None:
                if (bin_fields[0] in logs and logs[bin_fields[0]]) or (
                    bin_fields[1] in logs and logs[bin_fields[1]]
                ):
                    raise RuntimeError(
                        "CIC deposition is only implemented for linear-scaled axes"
                    )
            else:
                logs = {bin_fields[0]: False, bin_fields[1]: False}
            if any(accumulation) or fractional:
                raise RuntimeError(
                    "The accumulation and fractional keyword arguments must be "
                    "False for CIC deposition"
                )
        cls = ParticleProfile
    elif len(bin_fields) == 2:
        cls = Profile2D
    elif len(bin_fields) == 3:
        cls = Profile3D
    else:
        raise NotImplementedError
    if weight_field is not None and cls == ParticleProfile:
        (weight_field,) = data_source._determine_fields([weight_field])
        wf = data_source.ds._get_field_info(weight_field)
        if not wf.sampling_type == "particle":
            weight_field = None
    if not is_sequence(n_bins):
        n_bins = [n_bins] * len(bin_fields)
    if not is_sequence(accumulation):
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

    # Are the extrema all Nones? Then treat them as though extrema was set as None
    if extrema is None or not any(collapse(extrema.values())):
        ex = [
            data_source.quantities["Extrema"](f, non_zero=l)
            for f, l in zip(bin_fields, logs)
        ]
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
            except KeyError as e:
                try:
                    field_ex = list(extrema[bin_field])
                except KeyError:
                    raise RuntimeError(
                        "Could not find field {} or {} in extrema".format(
                            bin_field[-1], bin_field
                        )
                    ) from e

            if isinstance(field_ex[0], tuple):
                field_ex = [data_source.ds.quan(*f) for f in field_ex]
            if any([exi is None for exi in field_ex]):
                try:
                    ds_extrema = data_source.quantities.extrema(bin_field)
                except AttributeError:
                    # ytdata profile datasets don't have data_source.quantities
                    bf_vals = data_source[bin_field]
                    ds_extrema = data_source.ds.arr([bf_vals.min(), bf_vals.max()])
                for i, exi in enumerate(field_ex):
                    if exi is None:
                        field_ex[i] = ds_extrema[i]
                        # pad extrema by epsilon so cells at bin edges are
                        # not excluded
                        field_ex[i] -= (-1) ** i * np.spacing(field_ex[i])
            if units is not None and bin_field in units:
                for i, exi in enumerate(field_ex):
                    if hasattr(exi, "units"):
                        field_ex[i] = exi.to(units[bin_field])
                    else:
                        field_ex[i] = data_source.ds.quan(exi, units[bin_field])
                fe = data_source.ds.arr(field_ex)
            else:
                if hasattr(field_ex, "units"):
                    fe = field_ex.to(bf_units)
                else:
                    fe = data_source.ds.arr(field_ex, bf_units)
            fe.convert_to_units(bf_units)
            field_ex = [fe[0].v, fe[1].v]
            if is_sequence(field_ex[0]):
                field_ex[0] = data_source.ds.quan(field_ex[0][0], field_ex[0][1])
                field_ex[0] = field_ex[0].in_units(bf_units)
            if is_sequence(field_ex[1]):
                field_ex[1] = data_source.ds.quan(field_ex[1][0], field_ex[1][1])
                field_ex[1] = field_ex[1].in_units(bf_units)
            ex.append(field_ex)

    if override_bins is not None:
        o_bins = []
        for bin_field in bin_fields:
            bf_units = data_source.ds.field_info[bin_field].output_units
            try:
                field_obin = override_bins[bin_field[-1]]
            except KeyError:
                field_obin = override_bins[bin_field]

            if field_obin is None:
                o_bins.append(None)
                continue

            if isinstance(field_obin, tuple):
                field_obin = data_source.ds.arr(*field_obin)

            if units is not None and bin_field in units:
                fe = data_source.ds.arr(field_obin, units[bin_field])
            else:
                if hasattr(field_obin, "units"):
                    fe = field_obin.to(bf_units)
                else:
                    fe = data_source.ds.arr(field_obin, bf_units)
            fe.convert_to_units(bf_units)
            field_obin = fe.d
            o_bins.append(field_obin)

    args = [data_source]
    for f, n, (mi, ma), l in zip(bin_fields, n_bins, ex, logs):
        if mi <= 0 and l:
            raise YTIllDefinedBounds(mi, ma)
        args += [f, n, mi, ma, l]
    kwargs = dict(weight_field=weight_field)
    if cls is ParticleProfile:
        kwargs["deposition"] = deposition
    if override_bins is not None:
        for o_bin, ax in zip(o_bins, ["x", "y", "z"]):
            kwargs[f"override_bins_{ax}"] = o_bin
    obj = cls(*args, **kwargs)
    obj.accumulation = accumulation
    obj.fractional = fractional
    if fields is not None:
        obj.add_fields([field for field in fields])
    for field in fields:
        if fractional:
            obj.field_data[field] /= obj.field_data[field].sum()
        for axis, acc in enumerate(accumulation):
            if not acc:
                continue
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
                temp = (temp * temp_weight).cumsum(axis=0) / temp_weight.cumsum(axis=0)
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

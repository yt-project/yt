import abc
import weakref
from collections import defaultdict
from contextlib import contextmanager
from typing import List, Tuple, Union

import numpy as np

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.config import ytcfg
from yt.data_objects.field_data import YTFieldData
from yt.data_objects.profiles import create_profile
from yt.fields.field_exceptions import NeedsGridType
from yt.frontends.ytdata.utilities import save_as_dataset
from yt.funcs import get_output_filename, is_sequence, iter_fields, mylog
from yt.units._numpy_wrapper_functions import uconcatenate
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.utilities.exceptions import (
    YTCouldNotGenerateField,
    YTException,
    YTFieldNotFound,
    YTFieldNotParseable,
    YTFieldTypeNotFound,
    YTNonIndexedDataContainer,
    YTSpatialFieldUnitError,
)
from yt.utilities.object_registries import data_object_registry
from yt.utilities.on_demand_imports import _firefly as firefly
from yt.utilities.parameter_file_storage import ParameterFileStore


def sanitize_weight_field(ds, field, weight):
    field_object = ds._get_field_info(field)
    if weight is None:
        if field_object.sampling_type == "particle":
            if field_object.name[0] == "gas":
                ptype = ds._sph_ptypes[0]
            else:
                ptype = field_object.name[0]
            weight_field = (ptype, "particle_ones")
        else:
            weight_field = ("index", "ones")
    else:
        weight_field = weight
    return weight_field


def _get_ipython_key_completion(ds):
    # tuple-completion (ftype, fname) was added in IPython 8.0.0
    # with earlier versions, completion works with fname only
    # this implementation should work transparently with all IPython versions
    tuple_keys = ds.field_list + ds.derived_field_list
    fnames = list({k[1] for k in tuple_keys})
    return tuple_keys + fnames


class YTDataContainer(abc.ABC):
    """
    Generic YTDataContainer container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    and deal with passing back and forth field parameters.
    """

    _chunk_info = None
    _num_ghost_zones = 0
    _con_args: Tuple[str, ...] = ()
    _skip_add = False
    _container_fields: Tuple[Union[str, Tuple[str, str]], ...] = ()
    _tds_attrs: Tuple[str, ...] = ()
    _tds_fields: Tuple[str, ...] = ()
    _field_cache = None
    _index = None
    _key_fields: List[str]

    def __init__(self, ds, field_parameters):
        """
        Typically this is never called directly, but only due to inheritance.
        It associates a :class:`~yt.data_objects.static_output.Dataset` with the class,
        sets its initial set of fields, and the remainder of the arguments
        are passed as field_parameters.
        """
        # ds is typically set in the new object type created in
        # Dataset._add_object_class but it can also be passed as a parameter to the
        # constructor, in which case it will override the default.
        # This code ensures it is never not set.
        if ds is not None:
            self.ds = ds
        else:
            if not hasattr(self, "ds"):
                raise RuntimeError(
                    "Error: ds must be set either through class type "
                    "or parameter to the constructor"
                )

        self._current_particle_type = "all"
        self._current_fluid_type = self.ds.default_fluid_type
        self.ds.objects.append(weakref.proxy(self))
        mylog.debug("Appending object to %s (type: %s)", self.ds, type(self))
        self.field_data = YTFieldData()
        if self.ds.unit_system.has_current_mks:
            mag_unit = "T"
        else:
            mag_unit = "G"
        self._default_field_parameters = {
            "center": self.ds.arr(np.zeros(3, dtype="float64"), "cm"),
            "bulk_velocity": self.ds.arr(np.zeros(3, dtype="float64"), "cm/s"),
            "bulk_magnetic_field": self.ds.arr(np.zeros(3, dtype="float64"), mag_unit),
            "normal": self.ds.arr([0.0, 0.0, 1.0], ""),
        }
        if field_parameters is None:
            field_parameters = {}
        self._set_default_field_parameters()
        for key, val in field_parameters.items():
            self.set_field_parameter(key, val)

    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(*args, **kwargs)
        if hasattr(cls, "_type_name") and not cls._skip_add:
            name = getattr(cls, "_override_selector_name", cls._type_name)
            data_object_registry[name] = cls

    @property
    def pf(self):
        return getattr(self, "ds", None)

    @property
    def index(self):
        if self._index is not None:
            return self._index
        self._index = self.ds.index
        return self._index

    def _debug(self):
        """
        When called from within a derived field, this will run pdb.  However,
        during field detection, it will not.  This allows you to more easily
        debug fields that are being called on actual objects.
        """
        import pdb

        pdb.set_trace()

    def _set_default_field_parameters(self):
        self.field_parameters = {}
        for k, v in self._default_field_parameters.items():
            self.set_field_parameter(k, v)

    def _is_default_field_parameter(self, parameter):
        if parameter not in self._default_field_parameters:
            return False
        return (
            self._default_field_parameters[parameter]
            is self.field_parameters[parameter]
        )

    def apply_units(self, arr, units):
        try:
            arr.units.registry = self.ds.unit_registry
            return arr.to(units)
        except AttributeError:
            return self.ds.arr(arr, units=units)

    def _first_matching_field(self, field):
        for ftype, fname in self.ds.derived_field_list:
            if fname == field:
                return (ftype, fname)

        raise YTFieldNotFound(field, self.ds)

    def _set_center(self, center):
        if center is None:
            self.center = None
            return
        elif isinstance(center, YTArray):
            self.center = self.ds.arr(center.astype("float64"))
            self.center.convert_to_units("code_length")
        elif isinstance(center, (list, tuple, np.ndarray)):
            if isinstance(center[0], YTQuantity):
                self.center = self.ds.arr([c.copy() for c in center], dtype="float64")
                self.center.convert_to_units("code_length")
            else:
                self.center = self.ds.arr(center, "code_length", dtype="float64")
        elif isinstance(center, str):
            if center.lower() in ("c", "center"):
                self.center = self.ds.domain_center
            # is this dangerous for race conditions?
            elif center.lower() in ("max", "m"):
                self.center = self.ds.find_max(("gas", "density"))[1]
            elif center.startswith("max_"):
                field = self._first_matching_field(center[4:])
                self.center = self.ds.find_max(field)[1]
            elif center.lower() == "min":
                self.center = self.ds.find_min(("gas", "density"))[1]
            elif center.startswith("min_"):
                field = self._first_matching_field(center[4:])
                self.center = self.ds.find_min(field)[1]
        else:
            self.center = self.ds.arr(center, "code_length", dtype="float64")

        if self.center.ndim > 1:
            mylog.debug("Removing singleton dimensions from 'center'.")
            self.center = np.squeeze(self.center)
            if self.center.ndim > 1:
                msg = (
                    "center array must be 1 dimensional, supplied center has "
                    f"{self.center.ndim} dimensions with shape {self.center.shape}."
                )
                raise YTException(msg)

        self.set_field_parameter("center", self.center)

    def get_field_parameter(self, name, default=None):
        """
        This is typically only used by derived field functions, but
        it returns parameters used to generate fields.
        """
        if name in self.field_parameters:
            return self.field_parameters[name]
        else:
            return default

    def set_field_parameter(self, name, val):
        """
        Here we set up dictionaries that get passed up and down and ultimately
        to derived fields.
        """
        self.field_parameters[name] = val

    def has_field_parameter(self, name):
        """
        Checks if a field parameter is set.
        """
        return name in self.field_parameters

    def clear_data(self):
        """
        Clears out all data from the YTDataContainer instance, freeing memory.
        """
        self.field_data.clear()

    def has_key(self, key):
        """
        Checks if a data field already exists.
        """
        return key in self.field_data

    def keys(self):
        return self.field_data.keys()

    def _reshape_vals(self, arr):
        return arr

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        f = self._determine_fields([key])[0]
        if f not in self.field_data and key not in self.field_data:
            if f in self._container_fields:
                self.field_data[f] = self.ds.arr(self._generate_container_field(f))
                return self.field_data[f]
            else:
                self.get_data(f)
        # fi.units is the unit expression string. We depend on the registry
        # hanging off the dataset to define this unit object.
        # Note that this is less succinct so that we can account for the case
        # when there are, for example, no elements in the object.
        try:
            rv = self.field_data[f]
        except KeyError:
            if isinstance(f, tuple):
                fi = self.ds._get_field_info(*f)
            elif isinstance(f, bytes):
                fi = self.ds._get_field_info("unknown", f)
            rv = self.ds.arr(self.field_data[key], fi.units)
        return rv

    def _ipython_key_completions_(self):
        return _get_ipython_key_completion(self.ds)

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.field_data[key] = val

    def __delitem__(self, key):
        """
        Deletes a field
        """
        if key not in self.field_data:
            key = self._determine_fields(key)[0]
        del self.field_data[key]

    def _generate_field(self, field):
        ftype, fname = field
        finfo = self.ds._get_field_info(*field)
        with self._field_type_state(ftype, finfo):
            if fname in self._container_fields:
                tr = self._generate_container_field(field)
            if finfo.sampling_type == "particle":
                tr = self._generate_particle_field(field)
            else:
                tr = self._generate_fluid_field(field)
            if tr is None:
                raise YTCouldNotGenerateField(field, self.ds)
            return tr

    def _generate_fluid_field(self, field):
        # First we check the validator
        ftype, fname = field
        finfo = self.ds._get_field_info(ftype, fname)
        if self._current_chunk is None or self._current_chunk.chunk_type != "spatial":
            gen_obj = self
        else:
            gen_obj = self._current_chunk.objs[0]
            gen_obj.field_parameters = self.field_parameters
        try:
            finfo.check_available(gen_obj)
        except NeedsGridType as ngt_exception:
            rv = self._generate_spatial_fluid(field, ngt_exception.ghost_zones)
        else:
            rv = finfo(gen_obj)
        return rv

    def _generate_spatial_fluid(self, field, ngz):
        finfo = self.ds._get_field_info(*field)
        if finfo.units is None:
            raise YTSpatialFieldUnitError(field)
        units = finfo.units
        try:
            rv = self.ds.arr(np.zeros(self.ires.size, dtype="float64"), units)
            accumulate = False
        except YTNonIndexedDataContainer:
            # In this case, we'll generate many tiny arrays of unknown size and
            # then concatenate them.
            outputs = []
            accumulate = True
        ind = 0
        if ngz == 0:
            deps = self._identify_dependencies([field], spatial=True)
            deps = self._determine_fields(deps)
            for _io_chunk in self.chunks([], "io", cache=False):
                for _chunk in self.chunks([], "spatial", ngz=0, preload_fields=deps):
                    o = self._current_chunk.objs[0]
                    if accumulate:
                        rv = self.ds.arr(np.empty(o.ires.size, dtype="float64"), units)
                        outputs.append(rv)
                        ind = 0  # Does this work with mesh?
                    with o._activate_cache():
                        ind += o.select(
                            self.selector, source=self[field], dest=rv, offset=ind
                        )
        else:
            chunks = self.index._chunk(self, "spatial", ngz=ngz)
            for chunk in chunks:
                with self._chunked_read(chunk):
                    gz = self._current_chunk.objs[0]
                    gz.field_parameters = self.field_parameters
                    wogz = gz._base_grid
                    if accumulate:
                        rv = self.ds.arr(
                            np.empty(wogz.ires.size, dtype="float64"), units
                        )
                        outputs.append(rv)
                    ind += wogz.select(
                        self.selector,
                        source=gz[field][ngz:-ngz, ngz:-ngz, ngz:-ngz],
                        dest=rv,
                        offset=ind,
                    )
        if accumulate:
            rv = uconcatenate(outputs)
        return rv

    def _generate_particle_field(self, field):
        # First we check the validator
        ftype, fname = field
        if self._current_chunk is None or self._current_chunk.chunk_type != "spatial":
            gen_obj = self
        else:
            gen_obj = self._current_chunk.objs[0]
        try:
            finfo = self.ds._get_field_info(*field)
            finfo.check_available(gen_obj)
        except NeedsGridType as ngt_exception:
            if ngt_exception.ghost_zones != 0:
                raise NotImplementedError from ngt_exception
            size = self._count_particles(ftype)
            rv = self.ds.arr(np.empty(size, dtype="float64"), finfo.units)
            ind = 0
            for _io_chunk in self.chunks([], "io", cache=False):
                for _chunk in self.chunks(field, "spatial"):
                    x, y, z = (self[ftype, f"particle_position_{ax}"] for ax in "xyz")
                    if x.size == 0:
                        continue
                    mask = self._current_chunk.objs[0].select_particles(
                        self.selector, x, y, z
                    )
                    if mask is None:
                        continue
                    # This requests it from the grid and does NOT mask it
                    data = self[field][mask]
                    rv[ind : ind + data.size] = data
                    ind += data.size
        else:
            with self._field_type_state(ftype, finfo, gen_obj):
                rv = self.ds._get_field_info(*field)(gen_obj)
        return rv

    def _count_particles(self, ftype):
        for (f1, _f2), val in self.field_data.items():
            if f1 == ftype:
                return val.size
        size = 0
        for _io_chunk in self.chunks([], "io", cache=False):
            for _chunk in self.chunks([], "spatial"):
                x, y, z = (self[ftype, f"particle_position_{ax}"] for ax in "xyz")
                if x.size == 0:
                    continue
                size += self._current_chunk.objs[0].count_particles(
                    self.selector, x, y, z
                )
        return size

    def _generate_container_field(self, field):
        raise NotImplementedError

    def _parameter_iterate(self, seq):
        for obj in seq:
            old_fp = obj.field_parameters
            obj.field_parameters = self.field_parameters
            yield obj
            obj.field_parameters = old_fp

    def write_out(self, filename, fields=None, format="%0.16e"):
        """Write out the YTDataContainer object in a text file.

        This function will take a data object and produce a tab delimited text
        file containing the fields presently existing and the fields given in
        the ``fields`` list.

        Parameters
        ----------
        filename : String
            The name of the file to write to.

        fields : List of string, Default = None
            If this is supplied, these fields will be added to the list of
            fields to be saved to disk. If not supplied, whatever fields
            presently exist will be used.

        format : String, Default = "%0.16e"
            Format of numbers to be written in the file.

        Raises
        ------
        ValueError
            Raised when there is no existing field.

        YTException
            Raised when field_type of supplied fields is inconsistent with the
            field_type of existing fields.

        Examples
        --------
        >>> ds = fake_particle_ds()
        >>> sp = ds.sphere(ds.domain_center, 0.25)
        >>> sp.write_out("sphere_1.txt")
        >>> sp.write_out("sphere_2.txt", fields=["cell_volume"])
        """
        if fields is None:
            fields = sorted(self.field_data.keys())

        field_order = [("index", k) for k in self._key_fields]
        diff_fields = [field for field in fields if field not in field_order]
        field_order += diff_fields
        field_order = sorted(self._determine_fields(field_order))

        field_shapes = defaultdict(list)
        for field in field_order:
            shape = self[field].shape
            field_shapes[shape].append(field)

        # Check all fields have the same shape
        if len(field_shapes) != 1:
            err_msg = ["Got fields with different number of elements:\n"]
            for shape, these_fields in field_shapes.items():
                err_msg.append(f"\t {these_fields} with shape {shape}")
            raise YTException("\n".join(err_msg))

        with open(filename, "w") as fid:
            field_header = [str(f) for f in field_order]
            fid.write("\t".join(["#"] + field_header + ["\n"]))
            field_data = np.array([self.field_data[field] for field in field_order])
            for line in range(field_data.shape[1]):
                field_data[:, line].tofile(fid, sep="\t", format=format)
                fid.write("\n")

    def to_dataframe(self, fields):
        r"""Export a data object to a :class:`~pandas.DataFrame`.

        This function will take a data object and an optional list of fields
        and export them to a :class:`~pandas.DataFrame` object.
        If pandas is not importable, this will raise ImportError.

        Parameters
        ----------
        fields : list of strings or tuple field names
            This is the list of fields to be exported into
            the DataFrame.

        Returns
        -------
        df : :class:`~pandas.DataFrame`
            The data contained in the object.

        Examples
        --------
        >>> dd = ds.all_data()
        >>> df = dd.to_dataframe([("gas", "density"), ("gas", "temperature")])
        """
        from yt.utilities.on_demand_imports import _pandas as pd

        data = {}
        fields = self._determine_fields(fields)
        for field in fields:
            data[field[-1]] = self[field]
        df = pd.DataFrame(data)
        return df

    def to_astropy_table(self, fields):
        """
        Export region data to a :class:~astropy.table.table.QTable,
        which is a Table object which is unit-aware. The QTable can then
        be exported to an ASCII file, FITS file, etc.

        See the AstroPy Table docs for more details:
        http://docs.astropy.org/en/stable/table/

        Parameters
        ----------
        fields : list of strings or tuple field names
            This is the list of fields to be exported into
            the QTable.

        Examples
        --------
        >>> sp = ds.sphere("c", (1.0, "Mpc"))
        >>> t = sp.to_astropy_table([("gas", "density"), ("gas", "temperature")])
        """
        from astropy.table import QTable

        t = QTable()
        fields = self._determine_fields(fields)
        for field in fields:
            t[field[-1]] = self[field].to_astropy()
        return t

    def save_as_dataset(self, filename=None, fields=None):
        r"""Export a data object to a reloadable yt dataset.

        This function will take a data object and output a dataset
        containing either the fields presently existing or fields
        given in the ``fields`` list.  The resulting dataset can be
        reloaded as a yt dataset.

        Parameters
        ----------
        filename : str, optional
            The name of the file to be written.  If None, the name
            will be a combination of the original dataset and the type
            of data container.
        fields : list of string or tuple field names, optional
            If this is supplied, it is the list of fields to be saved to
            disk.  If not supplied, all the fields that have been queried
            will be saved.

        Returns
        -------
        filename : str
            The name of the file that has been created.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
        >>> sp = ds.sphere(ds.domain_center, (10, "Mpc"))
        >>> fn = sp.save_as_dataset(fields=[("gas", "density"), ("gas", "temperature")])
        >>> sphere_ds = yt.load(fn)
        >>> # the original data container is available as the data attribute
        >>> print(sds.data[("gas", "density")])
        [  4.46237613e-32   4.86830178e-32   4.46335118e-32 ...,   6.43956165e-30
           3.57339907e-30   2.83150720e-30] g/cm**3
        >>> ad = sphere_ds.all_data()
        >>> print(ad[("gas", "temperature")])
        [  1.00000000e+00   1.00000000e+00   1.00000000e+00 ...,   4.40108359e+04
           4.54380547e+04   4.72560117e+04] K

        """

        keyword = f"{str(self.ds)}_{self._type_name}"
        filename = get_output_filename(filename, keyword, ".h5")

        data = {}
        if fields is not None:
            for f in self._determine_fields(fields):
                data[f] = self[f]
        else:
            data.update(self.field_data)
        # get the extra fields needed to reconstruct the container
        tds_fields = tuple(("index", t) for t in self._tds_fields)
        for f in [f for f in self._container_fields + tds_fields if f not in data]:
            data[f] = self[f]
        data_fields = list(data.keys())

        need_grid_positions = False
        need_particle_positions = False
        ptypes = []
        ftypes = {}
        for field in data_fields:
            if field in self._container_fields:
                ftypes[field] = "grid"
                need_grid_positions = True
            elif self.ds.field_info[field].sampling_type == "particle":
                if field[0] not in ptypes:
                    ptypes.append(field[0])
                ftypes[field] = field[0]
                need_particle_positions = True
            else:
                ftypes[field] = "grid"
                need_grid_positions = True
        # projections and slices use px and py, so don't need positions
        if self._type_name in ["cutting", "proj", "slice", "quad_proj"]:
            need_grid_positions = False

        if need_particle_positions:
            for ax in self.ds.coordinates.axis_order:
                for ptype in ptypes:
                    p_field = (ptype, f"particle_position_{ax}")
                    if p_field in self.ds.field_info and p_field not in data:
                        data_fields.append(field)
                        ftypes[p_field] = p_field[0]
                        data[p_field] = self[p_field]
        if need_grid_positions:
            for ax in self.ds.coordinates.axis_order:
                g_field = ("index", ax)
                if g_field in self.ds.field_info and g_field not in data:
                    data_fields.append(g_field)
                    ftypes[g_field] = "grid"
                    data[g_field] = self[g_field]
                g_field = ("index", "d" + ax)
                if g_field in self.ds.field_info and g_field not in data:
                    data_fields.append(g_field)
                    ftypes[g_field] = "grid"
                    data[g_field] = self[g_field]

        extra_attrs = {
            arg: getattr(self, arg, None) for arg in self._con_args + self._tds_attrs
        }
        extra_attrs["con_args"] = repr(self._con_args)
        extra_attrs["data_type"] = "yt_data_container"
        extra_attrs["container_type"] = self._type_name
        extra_attrs["dimensionality"] = self._dimensionality
        save_as_dataset(
            self.ds, filename, data, field_types=ftypes, extra_attrs=extra_attrs
        )

        return filename

    def to_glue(self, fields, label="yt", data_collection=None):
        """
        Takes specific *fields* in the container and exports them to
        Glue (http://glueviz.org) for interactive
        analysis. Optionally add a *label*. If you are already within
        the Glue environment, you can pass a *data_collection* object,
        otherwise Glue will be started.
        """
        from glue.core import Data, DataCollection

        if ytcfg.get("yt", "internals", "within_testing"):
            from glue.core.application_base import Application as GlueApplication
        else:
            try:
                from glue.app.qt.application import GlueApplication
            except ImportError:
                from glue.qt.glue_application import GlueApplication
        gdata = Data(label=label)
        for component_name in fields:
            gdata.add_component(self[component_name], component_name)

        if data_collection is None:
            dc = DataCollection([gdata])
            app = GlueApplication(dc)
            try:
                app.start()
            except AttributeError:
                # In testing we're using a dummy glue application object
                # that doesn't have a start method
                pass
        else:
            data_collection.append(gdata)

    def create_firefly_object(
        self,
        datadir=None,
        fields_to_include=None,
        fields_units=None,
        default_decimation_factor=100,
        velocity_units="km/s",
        coordinate_units="kpc",
        show_unused_fields=0,
        *,
        JSONdir=None,
        **kwargs,
    ):
        r"""This function links a region of data stored in a yt dataset
        to the Python frontend API for [Firefly](http://github.com/ageller/Firefly),
        a browser-based particle visualization tool.

        Parameters
        ----------

        datadir : string
            Path to where any `.json` files should be saved. If a relative
            path will assume relative to `${HOME}`. A value of `None` will default to `${HOME}/Data`.

        fields_to_include : array_like of strings
            A list of fields that you want to include in your
            Firefly visualization for on-the-fly filtering and
            colormapping.

        default_decimation_factor : integer
            The factor by which you want to decimate each particle group
            by (e.g. if there are 1e7 total particles in your simulation
            you might want to set this to 100 at first). Randomly samples
            your data like `shuffled_data[::decimation_factor]` so as to
            not overtax a system. This is adjustable on a per particle group
            basis by changing the returned reader's
            `reader.particleGroup[i].decimation_factor` before calling
            `reader.writeToDisk()`.

        velocity_units : string
            The units that the velocity should be converted to in order to
            show streamlines in Firefly. Defaults to km/s.

        coordinate_units : string
            The units that the coordinates should be converted to. Defaults to
            kpc.

        show_unused_fields : boolean
            A flag to optionally print the fields that are available, in the
            dataset but were not explicitly requested to be tracked.

        Any additional keyword arguments are passed to
        firefly.data_reader.Reader.__init__

        Returns
        -------
        reader : Firefly.data_reader.Reader object
            A reader object from the Firefly, configured
            to output the current region selected

        Examples
        --------

            >>> ramses_ds = yt.load(
            ...     "/Users/agurvich/Desktop/yt_workshop/"
            ...     + "DICEGalaxyDisk_nonCosmological/output_00002/info_00002.txt"
            ... )

            >>> region = ramses_ds.sphere(ramses_ds.domain_center, (1000, "kpc"))

            >>> reader = region.create_firefly_object(
            ...     "IsoGalaxyRamses",
            ...     fields_to_include=[
            ...         "particle_extra_field_1",
            ...         "particle_extra_field_2",
            ...     ],
            ...     fields_units=["dimensionless", "dimensionless"],
            ... )

            >>> reader.settings["color"]["io"] = [1, 1, 0, 1]
            >>> reader.particleGroups[0].decimation_factor = 100
            >>> reader.writeToDisk()
        """

        ## handle default arguments
        if fields_to_include is None:
            fields_to_include = []
        if fields_units is None:
            fields_units = []

        ## handle input validation, if any
        if len(fields_units) != len(fields_to_include):
            raise RuntimeError("Each requested field must have units.")

        ## for safety, in case someone passes a float just cast it
        default_decimation_factor = int(default_decimation_factor)

        if JSONdir is not None:
            issue_deprecation_warning(
                "The 'JSONdir' keyword argument is a deprecated alias for 'datadir'."
                "Please use 'datadir' directly.",
                since="4.1",
            )
            datadir = JSONdir

        ## initialize a firefly reader instance
        reader = firefly.data_reader.Reader(
            datadir=datadir, clean_datadir=True, **kwargs
        )

        ## create a ParticleGroup object that contains *every* field
        for ptype in sorted(self.ds.particle_types_raw):
            ## skip this particle type if it has no particles in this dataset
            if self[ptype, "relative_particle_position"].shape[0] == 0:
                continue

            ## loop through the fields and print them to the screen
            if show_unused_fields:
                ## read the available extra fields from yt
                this_ptype_fields = self.ds.particle_fields_by_type[ptype]

                ## load the extra fields and print them
                for field in this_ptype_fields:
                    if field not in fields_to_include:
                        mylog.warning(
                            "detected (but did not request) %s %s", ptype, field
                        )

            field_arrays = []
            field_names = []

            ## explicitly go after the fields we want
            for field, units in zip(fields_to_include, fields_units):
                ## determine if you want to take the log of the field for Firefly
                log_flag = "log(" in units

                ## read the field array from the dataset
                this_field_array = self[ptype, field]

                ## fix the units string and prepend 'log' to the field for
                ##  the UI name
                if log_flag:
                    units = units[len("log(") : -1]
                    field = f"log{field}"

                ## perform the unit conversion and take the log if
                ##  necessary.
                this_field_array.in_units(units)
                if log_flag:
                    this_field_array = np.log10(this_field_array)

                ## add this array to the tracked arrays
                field_arrays += [this_field_array]
                field_names = np.append(field_names, [field], axis=0)

            ## flag whether we want to filter and/or color by these fields
            ##  we'll assume yes for both cases, this can be changed after
            ##  the reader object is returned to the user.
            field_filter_flags = np.ones(len(field_names))
            field_colormap_flags = np.ones(len(field_names))

            ## create a firefly ParticleGroup for this particle type
            pg = firefly.data_reader.ParticleGroup(
                UIname=ptype,
                coordinates=self[ptype, "relative_particle_position"].in_units(
                    coordinate_units
                ),
                velocities=self[ptype, "relative_particle_velocity"].in_units(
                    velocity_units
                ),
                field_arrays=field_arrays,
                field_names=field_names,
                field_filter_flags=field_filter_flags,
                field_colormap_flags=field_colormap_flags,
                decimation_factor=default_decimation_factor,
            )

            ## bind this particle group to the firefly reader object
            reader.addParticleGroup(pg)

        return reader

    # Numpy-like Operations
    def argmax(self, field, axis=None):
        r"""Return the values at which the field is maximized.

        This will, in a parallel-aware fashion, find the maximum value and then
        return to you the values at that maximum location that are requested
        for "axis".  By default it will return the spatial positions (in the
        natural coordinate system), but it can be any field

        Parameters
        ----------
        field : string or tuple field name
            The field to maximize.
        axis : string or list of strings, optional
            If supplied, the fields to sample along; if not supplied, defaults
            to the coordinate fields.  This can be the name of the coordinate
            fields (i.e., 'x', 'y', 'z') or a list of fields, but cannot be 0,
            1, 2.

        Returns
        -------
        A list of YTQuantities as specified by the axis argument.

        Examples
        --------

        >>> temp_at_max_rho = reg.argmax(
        ...     ("gas", "density"), axis=("gas", "temperature")
        ... )
        >>> max_rho_xyz = reg.argmax(("gas", "density"))
        >>> t_mrho, v_mrho = reg.argmax(
        ...     ("gas", "density"),
        ...     axis=[("gas", "temperature"), ("gas", "velocity_magnitude")],
        ... )
        >>> x, y, z = reg.argmax(("gas", "density"))

        """
        if axis is None:
            mv, pos0, pos1, pos2 = self.quantities.max_location(field)
            return pos0, pos1, pos2
        if isinstance(axis, str):
            axis = [axis]
        rv = self.quantities.sample_at_max_field_values(field, axis)
        if len(rv) == 2:
            return rv[1]
        return rv[1:]

    def argmin(self, field, axis=None):
        r"""Return the values at which the field is minimized.

        This will, in a parallel-aware fashion, find the minimum value and then
        return to you the values at that minimum location that are requested
        for "axis".  By default it will return the spatial positions (in the
        natural coordinate system), but it can be any field

        Parameters
        ----------
        field : string or tuple field name
            The field to minimize.
        axis : string or list of strings, optional
            If supplied, the fields to sample along; if not supplied, defaults
            to the coordinate fields.  This can be the name of the coordinate
            fields (i.e., 'x', 'y', 'z') or a list of fields, but cannot be 0,
            1, 2.

        Returns
        -------
        A list of YTQuantities as specified by the axis argument.

        Examples
        --------

        >>> temp_at_min_rho = reg.argmin(
        ...     ("gas", "density"), axis=("gas", "temperature")
        ... )
        >>> min_rho_xyz = reg.argmin(("gas", "density"))
        >>> t_mrho, v_mrho = reg.argmin(
        ...     ("gas", "density"),
        ...     axis=[("gas", "temperature"), ("gas", "velocity_magnitude")],
        ... )
        >>> x, y, z = reg.argmin(("gas", "density"))

        """
        if axis is None:
            mv, pos0, pos1, pos2 = self.quantities.min_location(field)
            return pos0, pos1, pos2
        if isinstance(axis, str):
            axis = [axis]
        rv = self.quantities.sample_at_min_field_values(field, axis)
        if len(rv) == 2:
            return rv[1]
        return rv[1:]

    def _compute_extrema(self, field):
        if self._extrema_cache is None:
            self._extrema_cache = {}
        if field not in self._extrema_cache:
            # Note we still need to call extrema for each field, as of right
            # now
            mi, ma = self.quantities.extrema(field)
            self._extrema_cache[field] = (mi, ma)
        return self._extrema_cache[field]

    _extrema_cache = None

    def max(self, field, axis=None):
        r"""Compute the maximum of a field, optionally along an axis.

        This will, in a parallel-aware fashion, compute the maximum of the
        given field.  Supplying an axis will result in a return value of a
        YTProjection, with method 'max' for maximum intensity.  If the max has
        already been requested, it will use the cached extrema value.

        Parameters
        ----------
        field : string or tuple field name
            The field to maximize.
        axis : string, optional
            If supplied, the axis to project the maximum along.

        Returns
        -------
        Either a scalar or a YTProjection.

        Examples
        --------

        >>> max_temp = reg.max(("gas", "temperature"))
        >>> max_temp_proj = reg.max(("gas", "temperature"), axis=("index", "x"))
        """
        if axis is None:
            rv = tuple(self._compute_extrema(f)[1] for f in iter_fields(field))
            if len(rv) == 1:
                return rv[0]
            return rv
        elif axis in self.ds.coordinates.axis_name:
            return self.ds.proj(field, axis, data_source=self, method="max")
        else:
            raise NotImplementedError(f"Unknown axis {axis}")

    def min(self, field, axis=None):
        r"""Compute the minimum of a field.

        This will, in a parallel-aware fashion, compute the minimum of the
        given field. Supplying an axis will result in a return value of a
        YTProjection, with method 'min' for minimum intensity.  If the min
        has already been requested, it will use the cached extrema value.

        Parameters
        ----------
        field : string or tuple field name
            The field to minimize.
        axis : string, optional
            If supplied, the axis to compute the minimum along.

        Returns
        -------
        Either a scalar or a YTProjection.

        Examples
        --------

        >>> min_temp = reg.min(("gas", "temperature"))
        >>> min_temp_proj = reg.min(("gas", "temperature"), axis=("index", "x"))
        """
        if axis is None:
            rv = tuple(self._compute_extrema(f)[0] for f in iter_fields(field))
            if len(rv) == 1:
                return rv[0]
            return rv
        elif axis in self.ds.coordinates.axis_name:
            return self.ds.proj(field, axis, data_source=self, method="min")
        else:
            raise NotImplementedError(f"Unknown axis {axis}")

    def std(self, field, axis=None, weight=None):
        """Compute the standard deviation of a field, optionally along
        an axis, with a weight.

        This will, in a parallel-ware fashion, compute the standard
        deviation of the given field. If an axis is supplied, it
        will return a projection, where the weight is also supplied.

        By default the weight field will be "ones" or "particle_ones",
        depending on the field, resulting in an unweighted standard
        deviation.

        Parameters
        ----------
        field : string or tuple field name
            The field to calculate the standard deviation of
        axis : string, optional
            If supplied, the axis to compute the standard deviation
            along (i.e., to project along)
        weight : string, optional
            The field to use as a weight.

        Returns
        -------
        Scalar or YTProjection.
        """
        weight_field = sanitize_weight_field(self.ds, field, weight)
        if axis in self.ds.coordinates.axis_name:
            r = self.ds.proj(
                field, axis, data_source=self, weight_field=weight_field, moment=2
            )
        elif axis is None:
            r = self.quantities.weighted_standard_deviation(field, weight_field)[0]
        else:
            raise NotImplementedError(f"Unknown axis {axis}")
        return r

    def ptp(self, field):
        r"""Compute the range of values (maximum - minimum) of a field.

        This will, in a parallel-aware fashion, compute the "peak-to-peak" of
        the given field.

        Parameters
        ----------
        field : string or tuple field name
            The field to average.

        Returns
        -------
        Scalar

        Examples
        --------

        >>> rho_range = reg.ptp(("gas", "density"))
        """
        ex = self._compute_extrema(field)
        return ex[1] - ex[0]

    def profile(
        self,
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
    ):
        r"""
        Create a 1, 2, or 3D profile object from this data_source.

        The dimensionality of the profile object is chosen by the number of
        fields given in the bin_fields argument.  This simply calls
        :func:`yt.data_objects.profiles.create_profile`.

        Parameters
        ----------
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
        >>> ad = ds.all_data()
        >>> profile = ad.profile(
        ...     ad,
        ...     [("gas", "density")],
        ...     [("gas", "temperature"), ("gas", "velocity_x")],
        ... )
        >>> print(profile.x)
        >>> print(profile["gas", "temperature"])
        >>> plot = profile.plot()
        """
        p = create_profile(
            self,
            bin_fields,
            fields,
            n_bins,
            extrema,
            logs,
            units,
            weight_field,
            accumulation,
            fractional,
            deposition,
        )
        return p

    def mean(self, field, axis=None, weight=None):
        r"""Compute the mean of a field, optionally along an axis, with a
        weight.

        This will, in a parallel-aware fashion, compute the mean of the
        given field.  If an axis is supplied, it will return a projection,
        where the weight is also supplied.  By default the weight field will be
        "ones" or "particle_ones", depending on the field being averaged,
        resulting in an unweighted average.

        Parameters
        ----------
        field : string or tuple field name
            The field to average.
        axis : string, optional
            If supplied, the axis to compute the mean along (i.e., to project
            along)
        weight : string, optional
            The field to use as a weight.

        Returns
        -------
        Scalar or YTProjection.

        Examples
        --------

        >>> avg_rho = reg.mean(("gas", "density"), weight="cell_volume")
        >>> rho_weighted_T = reg.mean(
        ...     ("gas", "temperature"), axis=("index", "y"), weight=("gas", "density")
        ... )
        """
        weight_field = sanitize_weight_field(self.ds, field, weight)
        if axis in self.ds.coordinates.axis_name:
            r = self.ds.proj(field, axis, data_source=self, weight_field=weight_field)
        elif axis is None:
            r = self.quantities.weighted_average_quantity(field, weight_field)
        else:
            raise NotImplementedError(f"Unknown axis {axis}")
        return r

    def sum(self, field, axis=None):
        r"""Compute the sum of a field, optionally along an axis.

        This will, in a parallel-aware fashion, compute the sum of the given
        field.  If an axis is specified, it will return a projection (using
        method type "sum", which does not take into account path length) along
        that axis.

        Parameters
        ----------
        field : string or tuple field name
            The field to sum.
        axis : string, optional
            If supplied, the axis to sum along.

        Returns
        -------
        Either a scalar or a YTProjection.

        Examples
        --------

        >>> total_vol = reg.sum("cell_volume")
        >>> cell_count = reg.sum(("index", "ones"), axis=("index", "x"))
        """
        # Because we're using ``sum`` to specifically mean a sum or a
        # projection with the method="sum", we do not utilize the ``mean``
        # function.
        if axis in self.ds.coordinates.axis_name:
            with self._field_parameter_state({"axis": axis}):
                r = self.ds.proj(field, axis, data_source=self, method="sum")
        elif axis is None:
            r = self.quantities.total_quantity(field)
        else:
            raise NotImplementedError(f"Unknown axis {axis}")
        return r

    def integrate(self, field, weight=None, axis=None, *, moment=1):
        r"""Compute the integral (projection) of a field along an axis.

        This projects a field along an axis.

        Parameters
        ----------
        field : string or tuple field name
            The field to project.
        weight : string or tuple field name
            The field to weight the projection by
        axis : string
            The axis to project along.
        moment : integer, optional
            for a weighted projection, moment = 1 (the default) corresponds to a
            weighted average. moment = 2 corresponds to a weighted standard
            deviation.

        Returns
        -------
        YTProjection

        Examples
        --------

        >>> column_density = reg.integrate(("gas", "density"), axis=("index", "z"))
        """
        if weight is not None:
            weight_field = sanitize_weight_field(self.ds, field, weight)
        else:
            weight_field = None
        if axis in self.ds.coordinates.axis_name:
            r = self.ds.proj(
                field, axis, data_source=self, weight_field=weight_field, moment=moment
            )
        else:
            raise NotImplementedError(f"Unknown axis {axis}")
        return r

    @property
    def _hash(self):
        s = f"{self}"
        try:
            import hashlib

            return hashlib.md5(s.encode("utf-8")).hexdigest()
        except ImportError:
            return s

    def __reduce__(self):
        args = tuple(
            [self.ds._hash(), self._type_name]
            + [getattr(self, n) for n in self._con_args]
            + [self.field_parameters]
        )
        return (_reconstruct_object, args)

    def clone(self):
        r"""Clone a data object.

        This will make a duplicate of a data object; note that the
        `field_parameters` may not necessarily be deeply-copied.  If you modify
        the field parameters in-place, it may or may not be shared between the
        objects, depending on the type of object that that particular field
        parameter is.

        Notes
        -----
        One use case for this is to have multiple identical data objects that
        are being chunked over in different orders.

        Examples
        --------

        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> sp = ds.sphere("c", 0.1)
        >>> sp_clone = sp.clone()
        >>> sp[("gas", "density")]
        >>> print(sp.field_data.keys())
        [("gas", "density")]
        >>> print(sp_clone.field_data.keys())
        []
        """
        args = self.__reduce__()
        return args[0](self.ds, *args[1][1:])

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = f"{self.__class__.__name__} ({self.ds}): "
        for i in self._con_args:
            try:
                s += ", {}={}".format(
                    i,
                    getattr(self, i).in_base(unit_system=self.ds.unit_system),
                )
            except AttributeError:
                s += f", {i}={getattr(self, i)}"
        return s

    @contextmanager
    def _field_parameter_state(self, field_parameters):
        # What we're doing here is making a copy of the incoming field
        # parameters, and then updating it with our own.  This means that we'll
        # be using our own center, if set, rather than the supplied one.  But
        # it also means that any additionally set values can override it.
        old_field_parameters = self.field_parameters
        new_field_parameters = field_parameters.copy()
        new_field_parameters.update(old_field_parameters)
        self.field_parameters = new_field_parameters
        yield
        self.field_parameters = old_field_parameters

    @contextmanager
    def _field_type_state(self, ftype, finfo, obj=None):
        if obj is None:
            obj = self
        old_particle_type = obj._current_particle_type
        old_fluid_type = obj._current_fluid_type
        fluid_types = self.ds.fluid_types
        if finfo.sampling_type == "particle" and ftype not in fluid_types:
            obj._current_particle_type = ftype
        else:
            obj._current_fluid_type = ftype
        yield
        obj._current_particle_type = old_particle_type
        obj._current_fluid_type = old_fluid_type

    def _tupleize_field(self, field):

        try:
            ftype, fname = field.name
            return ftype, fname
        except AttributeError:
            pass

        if is_sequence(field) and not isinstance(field, str):
            try:
                ftype, fname = field
                if not all(isinstance(_, str) for _ in field):
                    raise TypeError
                return ftype, fname
            except TypeError as e:
                raise YTFieldNotParseable(field) from e
            except ValueError:
                pass

        try:
            fname = field
            finfo = self.ds._get_field_info(field)
            if finfo.sampling_type == "particle":
                ftype = self._current_particle_type
                if hasattr(self.ds, "_sph_ptypes"):
                    ptypes = self.ds._sph_ptypes
                    if finfo.name[0] in ptypes:
                        ftype = finfo.name[0]
                    elif finfo.is_alias and finfo.alias_name[0] in ptypes:
                        ftype = self._current_fluid_type
            else:
                ftype = self._current_fluid_type
                if (ftype, fname) not in self.ds.field_info:
                    ftype = self.ds._last_freq[0]
            return ftype, fname
        except YTFieldNotFound:
            pass

        if isinstance(field, str):
            return "unknown", field

        raise YTFieldNotParseable(field)

    def _determine_fields(self, fields):
        if str(fields) in self.ds._determined_fields:
            return self.ds._determined_fields[str(fields)]
        explicit_fields = []
        for field in iter_fields(fields):
            if field in self._container_fields:
                explicit_fields.append(field)
                continue

            ftype, fname = self._tupleize_field(field)
            finfo = self.ds._get_field_info(ftype, fname)

            # really ugly check to ensure that this field really does exist somewhere,
            # in some naming convention, before returning it as a possible field type
            if (
                (ftype, fname) not in self.ds.field_info
                and (ftype, fname) not in self.ds.field_list
                and fname not in self.ds.field_list
                and (ftype, fname) not in self.ds.derived_field_list
                and fname not in self.ds.derived_field_list
                and (ftype, fname) not in self._container_fields
            ):
                raise YTFieldNotFound((ftype, fname), self.ds)

            # these tests are really insufficient as a field type may be valid, and the
            # field name may be valid, but not the combination (field type, field name)
            particle_field = finfo.sampling_type == "particle"
            local_field = finfo.local_sampling
            if local_field:
                pass
            elif particle_field and ftype not in self.ds.particle_types:
                raise YTFieldTypeNotFound(ftype, ds=self.ds)
            elif not particle_field and ftype not in self.ds.fluid_types:
                raise YTFieldTypeNotFound(ftype, ds=self.ds)
            explicit_fields.append((ftype, fname))

        self.ds._determined_fields[str(fields)] = explicit_fields
        return explicit_fields

    _tree = None

    @property
    def tiles(self):
        if self._tree is not None:
            return self._tree
        self._tree = AMRKDTree(self.ds, data_source=self)
        return self._tree

    @property
    def blocks(self):
        for _io_chunk in self.chunks([], "io"):
            for _chunk in self.chunks([], "spatial", ngz=0):
                # For grids this will be a grid object, and for octrees it will
                # be an OctreeSubset.  Note that we delegate to the sub-object.
                o = self._current_chunk.objs[0]
                cache_fp = o.field_parameters.copy()
                o.field_parameters.update(self.field_parameters)
                for b, m in o.select_blocks(self.selector):
                    if m is None:
                        continue
                    yield b, m
                o.field_parameters = cache_fp


# PR3124: Given that save_as_dataset is now the recommended method for saving
# objects (see Issue 2021 and references therein), the following has been re-written.
#
# Original comments (still true):
#
# In the future, this would be better off being set up to more directly
# reference objects or retain state, perhaps with a context manager.
#
# One final detail: time series or multiple datasets in a single pickle
# seems problematic.


def _get_ds_by_hash(hash):
    from yt.data_objects.static_output import Dataset

    if isinstance(hash, Dataset):
        return hash
    from yt.data_objects.static_output import _cached_datasets

    for ds in _cached_datasets.values():
        if ds._hash() == hash:
            return ds
    return None


def _reconstruct_object(*args, **kwargs):
    # returns a reconstructed YTDataContainer. As of PR 3124, we now return
    # the actual YTDataContainer rather than a (ds, YTDataContainer) tuple.

    # pull out some arguments
    dsid = args[0]  # the hash id
    dtype = args[1]  # DataContainer type (e.g., 'region')
    field_parameters = args[-1]  # the field parameters

    # re-instantiate the base dataset from the hash and ParameterFileStore
    ds = _get_ds_by_hash(dsid)
    override_weakref = False
    if not ds:
        override_weakref = True
        datasets = ParameterFileStore()
        ds = datasets.get_ds_hash(dsid)

    # instantiate the class with remainder of the args and adjust the state
    cls = getattr(ds, dtype)
    obj = cls(*args[2:-1])
    obj.field_parameters.update(field_parameters)

    # any nested ds references are weakref.proxy(ds), so need to ensure the ds
    # we just loaded persists when we leave this function (nosetests fail without
    # this) if we did not have an actual dataset as an argument.
    if hasattr(obj, "ds") and override_weakref:
        obj.ds = ds

    return obj

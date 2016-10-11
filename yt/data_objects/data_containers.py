"""
Various non-grid data containers.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import itertools
import uuid

import numpy as np
import weakref
import shelve

from collections import defaultdict
from contextlib import contextmanager

from yt.data_objects.particle_io import particle_handler_registry
from yt.fields.derived_field import \
    DerivedField
from yt.frontends.ytdata.utilities import \
    save_as_dataset
from yt.funcs import \
    get_output_filename, \
    mylog, \
    ensure_list, \
    fix_axis, \
    iterable
from yt.units.unit_object import UnitParseError
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
import yt.units.dimensions as ytdims
from yt.utilities.exceptions import \
    YTUnitConversionError, \
    YTFieldUnitError, \
    YTFieldUnitParseError, \
    YTSpatialFieldUnitError, \
    YTCouldNotGenerateField, \
    YTFieldNotParseable, \
    YTFieldNotFound, \
    YTFieldTypeNotFound, \
    YTDataSelectorNotImplemented, \
    YTDimensionalityError
from yt.utilities.lib.marching_cubes import \
    march_cubes_grid, march_cubes_grid_flux
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.parameter_file_storage import \
    ParameterFileStore
from yt.utilities.amr_kdtree.api import \
    AMRKDTree
from .derived_quantities import DerivedQuantityCollection
from yt.fields.field_exceptions import \
    NeedsGridType
import yt.geometry.selection_routines
from yt.geometry.selection_routines import \
    compose_selector
from yt.extern.six import add_metaclass, string_types

data_object_registry = {}

def force_array(item, shape):
    try:
        return item.copy()
    except AttributeError:
        if item:
            return np.ones(shape, dtype='bool')
        else:
            return np.zeros(shape, dtype='bool')

def restore_field_information_state(func):
    """
    A decorator that takes a function with the API of (self, grid, field)
    and ensures that after the function is called, the field_parameters will
    be returned to normal.
    """
    def save_state(self, grid, field=None, *args, **kwargs):
        old_params = grid.field_parameters
        grid.field_parameters = self.field_parameters
        tr = func(self, grid, field, *args, **kwargs)
        grid.field_parameters = old_params
        return tr
    return save_state

class YTFieldData(dict):
    """
    A Container object for field data, instead of just having it be a dict.
    """
    pass

class RegisteredDataContainer(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_type_name") and not cls._skip_add:
            data_object_registry[cls._type_name] = cls

@add_metaclass(RegisteredDataContainer)
class YTDataContainer(object):
    """
    Generic YTDataContainer container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    and deal with passing back and forth field parameters.
    """
    _chunk_info = None
    _num_ghost_zones = 0
    _con_args = ()
    _skip_add = False
    _container_fields = ()
    _tds_attrs = ()
    _tds_fields = ()
    _field_cache = None
    _index = None

    def __init__(self, ds, field_parameters):
        """
        Typically this is never called directly, but only due to inheritance.
        It associates a :class:`~yt.data_objects.api.Dataset` with the class,
        sets its initial set of fields, and the remainder of the arguments
        are passed as field_parameters.
        """
        # ds is typically set in the new object type created in Dataset._add_object_class
        # but it can also be passed as a parameter to the constructor, in which case it will
        # override the default. This code ensures it is never not set.
        if ds is not None:
            self.ds = ds
        else:
            if not hasattr(self, "ds"):
                raise RuntimeError("Error: ds must be set either through class type or parameter to the constructor")

        self._current_particle_type = "all"
        self._current_fluid_type = self.ds.default_fluid_type
        self.ds.objects.append(weakref.proxy(self))
        mylog.debug("Appending object to %s (type: %s)", self.ds, type(self))
        self.field_data = YTFieldData()
        self._default_field_parameters = {
            'center': self.ds.arr(np.zeros(3, dtype='float64'), 'cm'),
            'bulk_velocity': self.ds.arr(np.zeros(3, dtype='float64'), 'cm/s'),
            'normal': self.ds.arr([0.0, 0.0, 1.0], ''),
        }
        if field_parameters is None: field_parameters = {}
        self._set_default_field_parameters()
        for key, val in field_parameters.items():
            self.set_field_parameter(key, val)

    @property
    def pf(self):
        return getattr(self, 'ds', None)

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
        for k,v in self._default_field_parameters.items():
            self.set_field_parameter(k,v)

    def _is_default_field_parameter(self, parameter):
        if parameter not in self._default_field_parameters:
            return False
        return self._default_field_parameters[parameter] is \
          self.field_parameters[parameter]

    def apply_units(self, arr, units):
        return self.ds.arr(arr, input_units = units)

    def _set_center(self, center):
        if center is None:
            self.center = None
            return
        elif isinstance(center, YTArray):
            self.center = self.ds.arr(center.copy())
            self.center.convert_to_units('code_length')
        elif isinstance(center, (list, tuple, np.ndarray)):
            if isinstance(center[0], YTQuantity):
                self.center = self.ds.arr([c.copy() for c in center])
                self.center.convert_to_units('code_length')
            else:
                self.center = self.ds.arr(center, 'code_length')
        elif isinstance(center, string_types):
            if center.lower() in ("c", "center"):
                self.center = self.ds.domain_center
            # is this dangerous for race conditions?
            elif center.lower() in ("max", "m"):
                self.center = self.ds.find_max(("gas", "density"))[1]
            elif center.startswith("max_"):
                self.center = self.ds.find_max(center[4:])[1]
        else:
            self.center = self.ds.arr(center, 'code_length', dtype='float64')
        self.set_field_parameter('center', self.center)

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

    def convert(self, datatype):
        """
        This will attempt to convert a given unit to cgs from code units.
        It either returns the multiplicative factor or throws a KeyError.
        """
        return self.ds[datatype]

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
                self.field_data[f] = \
                    self.ds.arr(self._generate_container_field(f))
                return self.field_data[f]
            else:
                self.get_data(f)
        # fi.units is the unit expression string. We depend on the registry
        # hanging off the dataset to define this unit object.
        # Note that this is less succinct so that we can account for the case
        # when there are, for example, no elements in the object.
        rv = self.field_data.get(f, None)
        if rv is None:
            if isinstance(f, tuple):
                fi = self.ds._get_field_info(*f)
            elif isinstance(f, bytes):
                fi = self.ds._get_field_info("unknown", f)
            rv = self.ds.arr(self.field_data[key], fi.units)
        return rv

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
            if finfo.particle_type:
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
        if self._current_chunk is None or \
           self._current_chunk.chunk_type != "spatial":
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
        rv = self.ds.arr(np.empty(self.ires.size, dtype="float64"), units)
        ind = 0
        if ngz == 0:
            deps = self._identify_dependencies([field], spatial = True)
            deps = self._determine_fields(deps)
            for io_chunk in self.chunks([], "io", cache = False):
                for i,chunk in enumerate(self.chunks([], "spatial", ngz = 0,
                                                    preload_fields = deps)):
                    o = self._current_chunk.objs[0]
                    with o._activate_cache():
                        ind += o.select(self.selector, self[field], rv, ind)
        else:
            chunks = self.index._chunk(self, "spatial", ngz = ngz)
            for i, chunk in enumerate(chunks):
                with self._chunked_read(chunk):
                    gz = self._current_chunk.objs[0]
                    gz.field_parameters = self.field_parameters
                    wogz = gz._base_grid
                    ind += wogz.select(
                        self.selector,
                        gz[field][ngz:-ngz, ngz:-ngz, ngz:-ngz],
                        rv, ind)
        return rv

    def _generate_particle_field(self, field):
        # First we check the validator
        ftype, fname = field
        if self._current_chunk is None or \
           self._current_chunk.chunk_type != "spatial":
            gen_obj = self
        else:
            gen_obj = self._current_chunk.objs[0]
        try:
            finfo = self.ds._get_field_info(*field)
            finfo.check_available(gen_obj)
        except NeedsGridType as ngt_exception:
            if ngt_exception.ghost_zones != 0:
                raise NotImplementedError
            size = self._count_particles(ftype)
            rv = self.ds.arr(np.empty(size, dtype="float64"), finfo.units)
            ind = 0
            for io_chunk in self.chunks([], "io", cache = False):
                for i, chunk in enumerate(self.chunks(field, "spatial")):
                    x, y, z = (self[ftype, 'particle_position_%s' % ax]
                               for ax in 'xyz')
                    if x.size == 0: continue
                    mask = self._current_chunk.objs[0].select_particles(
                        self.selector, x, y, z)
                    if mask is None: continue
                    # This requests it from the grid and does NOT mask it
                    data = self[field][mask]
                    rv[ind:ind+data.size] = data
                    ind += data.size
        else:
            with self._field_type_state(ftype, finfo, gen_obj):
                rv = self.ds._get_field_info(*field)(gen_obj)
        return rv

    def _count_particles(self, ftype):
        for (f1, f2), val in self.field_data.items():
            if f1 == ftype:
                return val.size
        size = 0
        for io_chunk in self.chunks([], "io", cache = False):
            for i,chunk in enumerate(self.chunks([], "spatial")):
                x, y, z = (self[ftype, 'particle_position_%s' % ax]
                            for ax in 'xyz')
                if x.size == 0: continue
                size += self._current_chunk.objs[0].count_particles(
                    self.selector, x, y, z)
        return size

    def _generate_container_field(self, field):
        raise NotImplementedError

    def _parameter_iterate(self, seq):
        for obj in seq:
            old_fp = obj.field_parameters
            obj.field_parameters = self.field_parameters
            yield obj
            obj.field_parameters = old_fp

    _key_fields = None
    def write_out(self, filename, fields=None, format="%0.16e"):
        if fields is None: fields=sorted(self.field_data.keys())
        if self._key_fields is None: raise ValueError
        field_order = self._key_fields[:]
        for field in field_order: self[field]
        field_order += [field for field in fields if field not in field_order]
        fid = open(filename,"w")
        fid.write("\t".join(["#"] + field_order + ["\n"]))
        field_data = np.array([self.field_data[field] for field in field_order])
        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

    def save_object(self, name, filename = None):
        """
        Save an object.  If *filename* is supplied, it will be stored in
        a :mod:`shelve` file of that name.  Otherwise, it will be stored via
        :meth:`yt.data_objects.api.GridIndex.save_object`.
        """
        if filename is not None:
            ds = shelve.open(filename, protocol=-1)
            if name in ds:
                mylog.info("Overwriting %s in %s", name, filename)
            ds[name] = self
            ds.close()
        else:
            self.index.save_object(self, name)

    def to_dataframe(self, fields = None):
        r"""Export a data object to a pandas DataFrame.

        This function will take a data object and construct from it and
        optionally a list of fields a pandas DataFrame object.  If pandas is
        not importable, this will raise ImportError.

        Parameters
        ----------
        fields : list of strings or tuples, default None
            If this is supplied, it is the list of fields to be exported into
            the data frame.  If not supplied, whatever fields presently exist
            will be used.

        Returns
        -------
        df : DataFrame
            The data contained in the object.

        Examples
        --------

        >>> dd = ds.all_data()
        >>> df1 = dd.to_dataframe(["density", "temperature"])
        >>> dd["velocity_magnitude"]
        >>> df2 = dd.to_dataframe()
        """
        import pandas as pd
        data = {}
        if fields is not None:
            for f in fields:
                data[f] = self[f]
        else:
            data.update(self.field_data)
        df = pd.DataFrame(data)
        return df

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
        fields : list of strings or tuples, optional
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
        >>> fn = sp.save_as_dataset(fields=["density", "temperature"])
        >>> sphere_ds = yt.load(fn)
        >>> # the original data container is available as the data attribute
        >>> print (sds.data["density"])
        [  4.46237613e-32   4.86830178e-32   4.46335118e-32 ...,   6.43956165e-30
           3.57339907e-30   2.83150720e-30] g/cm**3
        >>> ad = sphere_ds.all_data()
        >>> print (ad["temperature"])
        [  1.00000000e+00   1.00000000e+00   1.00000000e+00 ...,   4.40108359e+04
           4.54380547e+04   4.72560117e+04] K

        """

        keyword = "%s_%s" % (str(self.ds), self._type_name)
        filename = get_output_filename(filename, keyword, ".h5")

        data = {}
        if fields is not None:
            for f in self._determine_fields(fields):
                data[f] = self[f]
        else:
            data.update(self.field_data)
        # get the extra fields needed to reconstruct the container
        tds_fields = tuple(self._determine_fields(list(self._tds_fields)))
        for f in [f for f in self._container_fields + tds_fields \
                  if f not in data]:
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
            elif self.ds.field_info[field].particle_type:
                if field[0] not in ptypes:
                    ptypes.append(field[0])
                ftypes[field] = field[0]
                need_particle_positions = True
            else:
                ftypes[field] = "grid"
                need_grid_positions = True
        # projections and slices use px and py, so don't need positions
        if self._type_name in ["cutting", "proj", "slice"]:
            need_grid_positions = False

        if need_particle_positions:
            for ax in "xyz":
                for ptype in ptypes:
                    p_field = (ptype, "particle_position_%s" % ax)
                    if p_field in self.ds.field_info and p_field not in data:
                        data_fields.append(field)
                        ftypes[p_field] = p_field[0]
                        data[p_field] = self[p_field]
        if need_grid_positions:
            for ax in "xyz":
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

        extra_attrs = dict([(arg, getattr(self, arg, None))
                            for arg in self._con_args + self._tds_attrs])
        extra_attrs["con_args"] = self._con_args
        extra_attrs["data_type"] = "yt_data_container"
        extra_attrs["container_type"] = self._type_name
        extra_attrs["dimensionality"] = self._dimensionality
        save_as_dataset(self.ds, filename, data, field_types=ftypes,
                        extra_attrs=extra_attrs)

        return filename

    def to_glue(self, fields, label="yt", data_collection=None):
        """
        Takes specific *fields* in the container and exports them to
        Glue (http://www.glueviz.org) for interactive
        analysis. Optionally add a *label*. If you are already within
        the Glue environment, you can pass a *data_collection* object,
        otherwise Glue will be started.
        """
        from glue.core import DataCollection, Data
        from glue.qt.glue_application import GlueApplication

        gdata = Data(label=label)
        for component_name in fields:
            gdata.add_component(self[component_name], component_name)

        if data_collection is None:
            dc = DataCollection([gdata])
            app = GlueApplication(dc)
            app.start()
        else:
            data_collection.append(gdata)

    # Numpy-like Operations
    def argmax(self, field, axis=None):
        r"""Return the values at which the field is maximized.

        This will, in a parallel-aware fashion, find the maximum value and then
        return to you the values at that maximum location that are requested
        for "axis".  By default it will return the spatial positions (in the
        natural coordinate system), but it can be any field

        Parameters
        ----------
        field : string or tuple of strings
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

        >>> temp_at_max_rho = reg.argmax("density", axis="temperature")
        >>> max_rho_xyz = reg.argmax("density")
        >>> t_mrho, v_mrho = reg.argmax("density", axis=["temperature",
        ...                 "velocity_magnitude"])
        >>> x, y, z = reg.argmax("density")

        """
        if axis is None:
            mv, pos0, pos1, pos2 = self.quantities.max_location(field)
            return pos0, pos1, pos2
        if isinstance(axis, string_types):
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
        field : string or tuple of strings
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

        >>> temp_at_min_rho = reg.argmin("density", axis="temperature")
        >>> min_rho_xyz = reg.argmin("density")
        >>> t_mrho, v_mrho = reg.argmin("density", axis=["temperature",
        ...                 "velocity_magnitude"])
        >>> x, y, z = reg.argmin("density")

        """
        if axis is None:
            mv, pos0, pos1, pos2 = self.quantities.min_location(field)
            return pos0, pos1, pos2
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
        YTProjection, with method 'mip' for maximum intensity.  If the max has
        already been requested, it will use the cached extrema value.

        Parameters
        ----------
        field : string or tuple of strings
            The field to maximize.
        axis : string, optional
            If supplied, the axis to project the maximum along.

        Returns
        -------
        Either a scalar or a YTProjection.

        Examples
        --------

        >>> max_temp = reg.max("temperature")
        >>> max_temp_proj = reg.max("temperature", axis="x")
        """
        if axis is None:
            rv = ()
            fields = ensure_list(field)
            for f in fields:
                rv += (self._compute_extrema(f)[1],)
            if len(fields) == 1:
                return rv[0]
            else:
                return rv
        elif axis in self.ds.coordinates.axis_name:
            r = self.ds.proj(field, axis, data_source=self, method="mip")
            return r
        else:
            raise NotImplementedError("Unknown axis %s" % axis)

    def min(self, field, axis=None):
        r"""Compute the minimum of a field.

        This will, in a parallel-aware fashion, compute the minimum of the
        given field.  Supplying an axis is not currently supported.  If the max
        has already been requested, it will use the cached extrema value.

        Parameters
        ----------
        field : string or tuple of strings
            The field to minimize.
        axis : string, optional
            If supplied, the axis to compute the minimum along.

        Returns
        -------
        Scalar.

        Examples
        --------

        >>> min_temp = reg.min("temperature")
        """
        if axis is None:
            rv = ()
            fields = ensure_list(field)
            for f in ensure_list(fields):
                rv += (self._compute_extrema(f)[0],)
            if len(fields) == 1:
                return rv[0]
            else:
                return rv
            return rv
        elif axis in self.ds.coordinates.axis_name:
            raise NotImplementedError("Minimum intensity projection not"
                                      " implemented.")
        else:
            raise NotImplementedError("Unknown axis %s" % axis)

    def std(self, field, weight=None):
        raise NotImplementedError

    def ptp(self, field):
        r"""Compute the range of values (maximum - minimum) of a field.

        This will, in a parallel-aware fashion, compute the "peak-to-peak" of
        the given field.

        Parameters
        ----------
        field : string or tuple of strings
            The field to average.

        Returns
        -------
        Scalar

        Examples
        --------

        >>> rho_range = reg.ptp("density")
        """
        ex = self._compute_extrema(field)
        return ex[1] - ex[0]

    def hist(self, field, weight = None, bins = None):
        raise NotImplementedError

    def mean(self, field, axis=None, weight='ones'):
        r"""Compute the mean of a field, optionally along an axis, with a
        weight.

        This will, in a parallel-aware fashion, compute the mean of the
        given field.  If an axis is supplied, it will return a projection,
        where the weight is also supplied.  By default the weight is "ones",
        resulting in a strict average.

        Parameters
        ----------
        field : string or tuple of strings
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

        >>> avg_rho = reg.mean("density", weight="cell_volume")
        >>> rho_weighted_T = reg.mean("temperature", axis="y", weight="density")
        """
        if axis in self.ds.coordinates.axis_name:
            r = self.ds.proj(field, axis, data_source=self, weight_field=weight)
        elif axis is None:
            if weight is None:
                r = self.quantities.total_quantity(field)
            else:
                r = self.quantities.weighted_average_quantity(field, weight)
        else:
            raise NotImplementedError("Unknown axis %s" % axis)
        return r

    def sum(self, field, axis=None):
        r"""Compute the sum of a field, optionally along an axis.

        This will, in a parallel-aware fashion, compute the sum of the given
        field.  If an axis is specified, it will return a projection (using
        method type "sum", which does not take into account path length) along
        that axis.

        Parameters
        ----------
        field : string or tuple of strings
            The field to sum.
        axis : string, optional
            If supplied, the axis to sum along.

        Returns
        -------
        Either a scalar or a YTProjection.

        Examples
        --------

        >>> total_vol = reg.sum("cell_volume")
        >>> cell_count = reg.sum("ones", axis="x")
        """
        # Because we're using ``sum`` to specifically mean a sum or a
        # projection with the method="sum", we do not utilize the ``mean``
        # function.
        if axis in self.ds.coordinates.axis_name:
            with self._field_parameter_state({'axis':axis}):
                r = self.ds.proj(field, axis, data_source=self, method="sum")
        elif axis is None:
            r = self.quantities.total_quantity(field)
        else:
            raise NotImplementedError("Unknown axis %s" % axis)
        return r

    def integrate(self, field, axis=None):
        r"""Compute the integral (projection) of a field along an axis.

        This projects a field along an axis.

        Parameters
        ----------
        field : string or tuple of strings
            The field to project.
        axis : string
            The axis to project along.

        Returns
        -------
        YTProjection

        Examples
        --------

        >>> column_density = reg.integrate("density", axis="z")
        """
        if axis in self.ds.coordinates.axis_name:
            r = self.ds.proj(field, axis, data_source=self)
        else:
            raise NotImplementedError("Unknown axis %s" % axis)
        return r

    @property
    def _hash(self):
        s = "%s" % self
        try:
            import hashlib
            return hashlib.md5(s.encode('utf-8')).hexdigest()
        except ImportError:
            return s

    def __reduce__(self):
        args = tuple([self.ds._hash(), self._type_name] +
                     [getattr(self, n) for n in self._con_args] +
                     [self.field_parameters])
        return (_reconstruct_object, args)

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = "%s (%s): " % (self.__class__.__name__, self.ds)
        for i in self._con_args:
            try:
                s += ", %s=%s" % (i, getattr(self, i).in_base(unit_system=self.ds.unit_system))
            except AttributeError:
                s += ", %s=%s" % (i, getattr(self, i))
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
    def _field_type_state(self, ftype, finfo, obj = None):
        if obj is None: obj = self
        old_particle_type = obj._current_particle_type
        old_fluid_type = obj._current_fluid_type
        if finfo.particle_type:
            obj._current_particle_type = ftype
        else:
            obj._current_fluid_type = ftype
        yield
        obj._current_particle_type = old_particle_type
        obj._current_fluid_type = old_fluid_type

    def _determine_fields(self, fields):
        fields = ensure_list(fields)
        explicit_fields = []
        for field in fields:
            if field in self._container_fields:
                explicit_fields.append(field)
                continue
            if isinstance(field, tuple):
                if len(field) != 2 or \
                   not isinstance(field[0], string_types) or \
                   not isinstance(field[1], string_types):
                    raise YTFieldNotParseable(field)
                ftype, fname = field
                finfo = self.ds._get_field_info(ftype, fname)
            elif isinstance(field, DerivedField):
                ftype, fname = field.name
                finfo = field
            else:
                fname = field
                finfo = self.ds._get_field_info("unknown", fname)
                if finfo.particle_type:
                    ftype = self._current_particle_type
                else:
                    ftype = self._current_fluid_type
                    if (ftype, fname) not in self.ds.field_info:
                        ftype = self.ds._last_freq[0]

                # really ugly check to ensure that this field really does exist somewhere,
                # in some naming convention, before returning it as a possible field type
                if (ftype,fname) not in self.ds.field_info and \
                        (ftype,fname) not in self.ds.field_list and \
                        fname not in self.ds.field_list and \
                        (ftype,fname) not in self.ds.derived_field_list and \
                        fname not in self.ds.derived_field_list and \
                        (ftype,fname) not in self._container_fields:
                    raise YTFieldNotFound((ftype,fname),self.ds)

            # these tests are really insufficient as a field type may be valid, and the
            # field name may be valid, but not the combination (field type, field name)
            if finfo.particle_type and ftype not in self.ds.particle_types:
                raise YTFieldTypeNotFound(ftype, ds=self.ds)
            elif not finfo.particle_type and ftype not in self.ds.fluid_types:
                raise YTFieldTypeNotFound(ftype, ds=self.ds)
            explicit_fields.append((ftype, fname))
        return explicit_fields

    _tree = None

    @property
    def tiles(self):
        if self._tree is not None: return self._tree
        self._tree = AMRKDTree(self.ds, data_source=self)
        return self._tree

    @property
    def blocks(self):
        for io_chunk in self.chunks([], "io"):
            for i,chunk in enumerate(self.chunks([], "spatial", ngz = 0)):
                # For grids this will be a grid object, and for octrees it will
                # be an OctreeSubset.  Note that we delegate to the sub-object.
                o = self._current_chunk.objs[0]
                for b, m in o.select_blocks(self.selector):
                    if m is None: continue
                    yield b, m

class GenerationInProgress(Exception):
    def __init__(self, fields):
        self.fields = fields
        super(GenerationInProgress, self).__init__()

class YTSelectionContainer(YTDataContainer, ParallelAnalysisInterface):
    _locked = False
    _sort_by = None
    _selector = None
    _current_chunk = None
    _data_source = None
    _dimensionality = None

    def __init__(self, ds, field_parameters, data_source=None):
        ParallelAnalysisInterface.__init__(self)
        super(YTSelectionContainer, self).__init__(ds, field_parameters)
        self._data_source = data_source
        if data_source is not None:
            if data_source.ds is not self.ds:
                raise RuntimeError("Attempted to construct a DataContainer with a data_source "
                                   "from a different DataSet", ds, data_source.ds)
            if data_source._dimensionality < self._dimensionality:
                raise RuntimeError("Attempted to construct a DataContainer with a data_source "
                                   "of lower dimensionality (%u vs %u)" %
                                    (data_source._dimensionality, self._dimensionality))
            self.field_parameters.update(data_source.field_parameters)
        self.quantities = DerivedQuantityCollection(self)

    @property
    def selector(self):
        if self._selector is not None: return self._selector
        s_module = getattr(self, '_selector_module',
                           yt.geometry.selection_routines)
        sclass = getattr(s_module,
                         "%s_selector" % self._type_name, None)
        if sclass is None:
            raise YTDataSelectorNotImplemented(self._type_name)

        if self._data_source is not None:
            self._selector = compose_selector(self, self._data_source.selector, sclass(self))
        else:
            self._selector = sclass(self)
        return self._selector

    def chunks(self, fields, chunking_style, **kwargs):
        # This is an iterator that will yield the necessary chunks.
        self.get_data() # Ensure we have built ourselves
        if fields is None: fields = []
        for chunk in self.index._chunk(self, chunking_style, **kwargs):
            with self._chunked_read(chunk):
                self.get_data(fields)
                # NOTE: we yield before releasing the context
                yield self

    def _identify_dependencies(self, fields_to_get, spatial = False):
        inspected = 0
        fields_to_get = fields_to_get[:]
        for field in itertools.cycle(fields_to_get):
            if inspected >= len(fields_to_get): break
            inspected += 1
            fi = self.ds._get_field_info(*field)
            fd = self.ds.field_dependencies.get(field, None) or \
                 self.ds.field_dependencies.get(field[1], None)
            # This is long overdue.  Any time we *can't* find a field
            # dependency -- for instance, if the derived field has been added
            # after dataset instantiation -- let's just try to
            # recalculate it.
            if fd is None:
                try:
                    fd = fi.get_dependencies(ds = self.ds)
                    self.ds.field_dependencies[field] = fd
                except:
                    continue
            requested = self._determine_fields(list(set(fd.requested)))
            deps = [d for d in requested if d not in fields_to_get]
            fields_to_get += deps
        return fields_to_get

    def get_data(self, fields=None):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if fields is None: return
        nfields = []
        apply_fields = defaultdict(list)
        for field in self._determine_fields(fields):
            if field[0] in self.ds.filtered_particle_types:
                f = self.ds.known_filters[field[0]]
                apply_fields[field[0]].append(
                    (f.filtered_type, field[1]))
            else:
                nfields.append(field)
        for filter_type in apply_fields:
            f = self.ds.known_filters[filter_type]
            with f.apply(self):
                self.get_data(apply_fields[filter_type])
        fields = nfields
        if len(fields) == 0: return
        # Now we collect all our fields
        # Here is where we need to perform a validation step, so that if we
        # have a field requested that we actually *can't* yet get, we put it
        # off until the end.  This prevents double-reading fields that will
        # need to be used in spatial fields later on.
        fields_to_get = []
        # This will be pre-populated with spatial fields
        fields_to_generate = []
        for field in self._determine_fields(fields):
            if field in self.field_data: continue
            finfo = self.ds._get_field_info(*field)
            try:
                finfo.check_available(self)
            except NeedsGridType:
                fields_to_generate.append(field)
                continue
            fields_to_get.append(field)
        if len(fields_to_get) == 0 and len(fields_to_generate) == 0:
            return
        elif self._locked is True:
            raise GenerationInProgress(fields)
        # Track which ones we want in the end
        ofields = set(list(self.field_data.keys())
                    + fields_to_get
                    + fields_to_generate)
        # At this point, we want to figure out *all* our dependencies.
        fields_to_get = self._identify_dependencies(fields_to_get,
            self._spatial)
        # We now split up into readers for the types of fields
        fluids, particles = [], []
        finfos = {}
        for ftype, fname in fields_to_get:
            finfo = self.ds._get_field_info(ftype, fname)
            finfos[ftype, fname] = finfo
            if finfo.particle_type:
                particles.append((ftype, fname))
            elif (ftype, fname) not in fluids:
                fluids.append((ftype, fname))
        # The _read method will figure out which fields it needs to get from
        # disk, and return a dict of those fields along with the fields that
        # need to be generated.
        read_fluids, gen_fluids = self.index._read_fluid_fields(
                                        fluids, self, self._current_chunk)
        for f, v in read_fluids.items():
            self.field_data[f] = self.ds.arr(v, input_units = finfos[f].units)
            self.field_data[f].convert_to_units(finfos[f].output_units)

        read_particles, gen_particles = self.index._read_particle_fields(
                                        particles, self, self._current_chunk)
        for f, v in read_particles.items():
            self.field_data[f] = self.ds.arr(v, input_units = finfos[f].units)
            self.field_data[f].convert_to_units(finfos[f].output_units)

        fields_to_generate += gen_fluids + gen_particles
        self._generate_fields(fields_to_generate)
        for field in list(self.field_data.keys()):
            if field not in ofields:
                self.field_data.pop(field)

    def _generate_fields(self, fields_to_generate):
        index = 0
        with self._field_lock():
            # At this point, we assume that any fields that are necessary to
            # *generate* a field are in fact already available to us.  Note
            # that we do not make any assumption about whether or not the
            # fields have a spatial requirement.  This will be checked inside
            # _generate_field, at which point additional dependencies may
            # actually be noted.
            while any(f not in self.field_data for f in fields_to_generate):
                field = fields_to_generate[index % len(fields_to_generate)]
                index += 1
                if field in self.field_data: continue
                fi = self.ds._get_field_info(*field)
                try:
                    fd = self._generate_field(field)
                    if fd is None:
                        raise RuntimeError
                    if fi.units is None:
                        # first time calling a field with units='auto', so we
                        # infer the units from the units of the data we get back
                        # from the field function and use these units for future
                        # field accesses
                        units = getattr(fd, 'units', '')
                        if units == '':
                            dimensions = ytdims.dimensionless
                        else:
                            dimensions = units.dimensions
                            units = str(units.get_base_equivalent(self.ds.unit_system.name))
                        if fi.dimensions != dimensions:
                            raise YTDimensionalityError(fi.dimensions, dimensions)
                        fi.units = units
                        self.field_data[field] = self.ds.arr(fd, units)
                        msg = ("Field %s was added without specifying units, "
                               "assuming units are %s")
                        mylog.warn(msg % (fi.name, units))
                    try:
                        fd.convert_to_units(fi.units)
                    except AttributeError:
                        # If the field returns an ndarray, coerce to a
                        # dimensionless YTArray and verify that field is
                        # supposed to be unitless
                        fd = self.ds.arr(fd, '')
                        if fi.units != '':
                            raise YTFieldUnitError(fi, fd.units)
                    except YTUnitConversionError:
                        raise YTFieldUnitError(fi, fd.units)
                    except UnitParseError:
                        raise YTFieldUnitParseError(fi)
                    self.field_data[field] = fd
                except GenerationInProgress as gip:
                    for f in gip.fields:
                        if f not in fields_to_generate:
                            fields_to_generate.append(f)

    @contextmanager
    def _field_lock(self):
        self._locked = True
        yield
        self._locked = False

    @contextmanager
    def _chunked_read(self, chunk):
        # There are several items that need to be swapped out
        # field_data, size, shape
        old_field_data, self.field_data = self.field_data, YTFieldData()
        old_chunk, self._current_chunk = self._current_chunk, chunk
        old_locked, self._locked = self._locked, False
        yield
        self.field_data = old_field_data
        self._current_chunk = old_chunk
        self._locked = old_locked

    @contextmanager
    def _activate_cache(self):
        cache = self._field_cache or {}
        old_fields = {}
        for field in (f for f in cache if f in self.field_data):
            old_fields[field] = self.field_data[field]
        self.field_data.update(cache)
        yield
        for field in cache:
            self.field_data.pop(field)
            if field in old_fields:
                self.field_data[field] = old_fields.pop(field)
        self._field_cache = None

    def _initialize_cache(self, cache):
        # Wipe out what came before
        self._field_cache = {}
        self._field_cache.update(cache)

    @property
    def icoords(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.icoords

    @property
    def fcoords(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.fcoords

    @property
    def ires(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.ires

    @property
    def fwidth(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.fwidth

    @property
    def fcoords_vertex(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.fcoords_vertex


class YTSelectionContainer0D(YTSelectionContainer):
    _spatial = False
    _dimensionality = 0
    def __init__(self, ds, field_parameters = None, data_source = None):
        super(YTSelectionContainer0D, self).__init__(
            ds, field_parameters, data_source)

class YTSelectionContainer1D(YTSelectionContainer):
    _spatial = False
    _dimensionality = 1
    def __init__(self, ds, field_parameters = None, data_source = None):
        super(YTSelectionContainer1D, self).__init__(
            ds, field_parameters, data_source)
        self._grids = None
        self._sortkey = None
        self._sorted = {}

class YTSelectionContainer2D(YTSelectionContainer):
    _key_fields = ['px','py','pdx','pdy']
    _dimensionality = 2
    """
    Prepares the YTSelectionContainer2D, normal to *axis*.  If *axis* is 4, we are not
    aligned with any axis.
    """
    _spatial = False
    def __init__(self, axis, ds, field_parameters = None, data_source = None):
        super(YTSelectionContainer2D, self).__init__(
            ds, field_parameters, data_source)
        # We need the ds, which will exist by now, for fix_axis.
        self.axis = fix_axis(axis, self.ds)
        self.set_field_parameter("axis", axis)

    def _convert_field_name(self, field):
        return field

    def _get_pw(self, fields, center, width, origin, plot_type):
        from yt.visualization.plot_window import \
            get_window_parameters, PWViewerMPL
        from yt.visualization.fixed_resolution import \
            FixedResolutionBuffer as frb
        axis = self.axis
        skip = self._key_fields
        skip += list(set(frb._exclude_fields).difference(set(self._key_fields)))
        self.fields = [k for k in self.field_data if k not in skip]
        if fields is not None:
            self.fields = ensure_list(fields) + self.fields
        if len(self.fields) == 0:
            raise ValueError("No fields found to plot in get_pw")
        (bounds, center, display_center) = \
            get_window_parameters(axis, center, width, self.ds)
        pw = PWViewerMPL(self, bounds, fields=self.fields, origin=origin,
                         frb_generator=frb, plot_type=plot_type)
        pw._setup_plots()
        return pw


    def to_frb(self, width, resolution, center=None, height=None,
               periodic = False):
        r"""This function returns a FixedResolutionBuffer generated from this
        object.

        A FixedResolutionBuffer is an object that accepts a variable-resolution
        2D object and transforms it into an NxM bitmap that can be plotted,
        examined or processed.  This is a convenience function to return an FRB
        directly from an existing 2D data object.

        Parameters
        ----------
        width : width specifier
            This can either be a floating point value, in the native domain
            units of the simulation, or a tuple of the (value, unit) style.
            This will be the width of the FRB.
        height : height specifier
            This will be the physical height of the FRB, by default it is equal
            to width.  Note that this will not make any corrections to
            resolution for the aspect ratio.
        resolution : int or tuple of ints
            The number of pixels on a side of the final FRB.  If iterable, this
            will be the width then the height.
        center : array-like of floats, optional
            The center of the FRB.  If not specified, defaults to the center of
            the current object.
        periodic : bool
            Should the returned Fixed Resolution Buffer be periodic?  (default:
            False).

        Returns
        -------
        frb : :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`
            A fixed resolution buffer, which can be queried for fields.

        Examples
        --------

        >>> proj = ds.proj("Density", 0)
        >>> frb = proj.to_frb( (100.0, 'kpc'), 1024)
        >>> write_image(np.log10(frb["Density"]), 'density_100kpc.png')
        """

        if (self.ds.geometry == "cylindrical" and self.axis == 1) or \
            (self.ds.geometry == "polar" and self.axis == 2):
            if center is not None and center != (0.0, 0.0):
                raise NotImplementedError(
                    "Currently we only support images centered at R=0. " +
                    "We plan to generalize this in the near future")
            from yt.visualization.fixed_resolution import CylindricalFixedResolutionBuffer
            if iterable(width):
                radius = max(width)
            else:
                radius = width
            if iterable(resolution): resolution = max(resolution)
            frb = CylindricalFixedResolutionBuffer(self, radius, resolution)
            return frb

        if center is None:
            center = self.center
            if center is None:
                center = (self.ds.domain_right_edge
                        + self.ds.domain_left_edge)/2.0
        elif iterable(center) and not isinstance(center, YTArray):
            center = self.ds.arr(center, 'code_length')
        if iterable(width):
            w, u = width
            if isinstance(w, tuple) and isinstance(u, tuple):
                height = u
                w, u = w
            width = self.ds.quan(w, input_units = u)
        elif not isinstance(width, YTArray):
            width = self.ds.quan(width, 'code_length')
        if height is None:
            height = width
        elif iterable(height):
            h, u = height
            height = self.ds.quan(h, input_units = u)
        if not iterable(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import FixedResolutionBuffer
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        bounds = (center[xax] - width*0.5, center[xax] + width*0.5,
                  center[yax] - height*0.5, center[yax] + height*0.5)
        frb = FixedResolutionBuffer(self, bounds, resolution,
                                    periodic = periodic)
        return frb

class YTSelectionContainer3D(YTSelectionContainer):
    """
    Returns an instance of YTSelectionContainer3D, or prepares one.  Usually only
    used as a base class.  Note that *center* is supplied, but only used
    for fields and quantities that require it.
    """
    _key_fields = ['x','y','z','dx','dy','dz']
    _spatial = False
    _num_ghost_zones = 0
    _dimensionality = 3
    def __init__(self, center, ds, field_parameters = None, data_source = None):
        super(YTSelectionContainer3D, self).__init__(ds, field_parameters, data_source)
        self._set_center(center)
        self.coords = None
        self._grids = None

    def cut_region(self, field_cuts, field_parameters=None):
        """
        Return a YTCutRegion, where the a cell is identified as being inside
        the cut region based on the value of one or more fields.  Note that in
        previous versions of yt the name 'grid' was used to represent the data
        object used to construct the field cut, as of yt 3.0, this has been
        changed to 'obj'.

        Parameters
        ----------
        field_cuts : list of strings
           A list of conditionals that will be evaluated. In the namespace
           available, these conditionals will have access to 'obj' which is a
           data object of unknown shape, and they must generate a boolean array.
           For instance, conditionals = ["obj['temperature'] < 1e3"]
        field_parameters : dictionary
           A dictionary of field parameters to be used when applying the field
           cuts.

        Examples
        --------
        To find the total mass of hot gas with temperature greater than 10^6 K
        in your volume:

        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.cut_region(["obj['temperature'] > 1e6"])
        >>> print cr.quantities.total_quantity("cell_mass").in_units('Msun')
        """
        cr = self.ds.cut_region(self, field_cuts,
                                field_parameters=field_parameters)
        return cr

    def extract_isocontours(self, field, value, filename = None,
                            rescale = False, sample_values = None):
        r"""This identifies isocontours on a cell-by-cell basis, with no
        consideration of global connectedness, and returns the vertices of the
        Triangles in that isocontour.

        This function simply returns the vertices of all the triangles
        calculated by the `marching cubes
        <http://en.wikipedia.org/wiki/Marching_cubes>`_ algorithm; for more
        complex operations, such as identifying connected sets of cells above a
        given threshold, see the extract_connected_sets function.  This is more
        useful for calculating, for instance, total isocontour area, or
        visualizing in an external program (such as `MeshLab
        <http://meshlab.sf.net>`_.)

        Parameters
        ----------
        field : string
            Any field that can be obtained in a data object.  This is the field
            which will be isocontoured.
        value : float
            The value at which the isocontour should be calculated.
        filename : string, optional
            If supplied, this file will be filled with the vertices in .obj
            format.  Suitable for loading into meshlab.
        rescale : bool, optional
            If true, the vertices will be rescaled within their min/max.
        sample_values : string, optional
            Any field whose value should be extracted at the center of each
            triangle.

        Returns
        -------
        verts : array of floats
            The array of vertices, x,y,z.  Taken in threes, these are the
            triangle vertices.
        samples : array of floats
            If `sample_values` is specified, this will be returned and will
            contain the values of the field specified at the center of each
            triangle.

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        output the vertices to "triangles.obj" after rescaling them.

        >>> dd = ds.all_data()
        >>> rho = dd.quantities["WeightedAverageQuantity"](
        ...     "Density", weight="CellMassMsun")
        >>> verts = dd.extract_isocontours("Density", rho,
        ...             "triangles.obj", True)
        """
        verts = []
        samples = []
        for block, mask in self.blocks:
            my_verts = self._extract_isocontours_from_grid(
                block, mask, field, value, sample_values)
            if sample_values is not None:
                my_verts, svals = my_verts
                samples.append(svals)
            verts.append(my_verts)
        verts = np.concatenate(verts).transpose()
        verts = self.comm.par_combine_object(verts, op='cat', datatype='array')
        verts = verts.transpose()
        if sample_values is not None:
            samples = np.concatenate(samples)
            samples = self.comm.par_combine_object(samples, op='cat',
                                datatype='array')
        if rescale:
            mi = np.min(verts, axis=0)
            ma = np.max(verts, axis=0)
            verts = (verts - mi) / (ma - mi).max()
        if filename is not None and self.comm.rank == 0:
            if hasattr(filename, "write"): f = filename
            else: f = open(filename, "w")
            for v1 in verts:
                f.write("v %0.16e %0.16e %0.16e\n" % (v1[0], v1[1], v1[2]))
            for i in range(len(verts)/3):
                f.write("f %s %s %s\n" % (i*3+1, i*3+2, i*3+3))
            if not hasattr(filename, "write"): f.close()
        if sample_values is not None:
            return verts, samples
        return verts

    def _extract_isocontours_from_grid(self, grid, mask, field, value,
                                       sample_values=None):
        vc_fields = [field]
        if sample_values is not None:
            vc_fields.append(sample_values)

        vc_data = grid.get_vertex_centered_data(vc_fields, no_ghost=False)
        try:
            svals = vc_data[sample_values]
        except KeyError:
            svals = None

        my_verts = march_cubes_grid(value, vc_data[field], mask,
            grid.LeftEdge, grid.dds, svals)
        return my_verts

    def calculate_isocontour_flux(self, field, value,
                    field_x, field_y, field_z, fluxing_field = None):
        r"""This identifies isocontours on a cell-by-cell basis, with no
        consideration of global connectedness, and calculates the flux over
        those contours.

        This function will conduct `marching cubes
        <http://en.wikipedia.org/wiki/Marching_cubes>`_ on all the cells in a
        given data container (grid-by-grid), and then for each identified
        triangular segment of an isocontour in a given cell, calculate the
        gradient (i.e., normal) in the isocontoured field, interpolate the local
        value of the "fluxing" field, the area of the triangle, and then return:

        area * local_flux_value * (n dot v)

        Where area, local_value, and the vector v are interpolated at the barycenter
        (weighted by the vertex values) of the triangle.  Note that this
        specifically allows for the field fluxing across the surface to be
        *different* from the field being contoured.  If the fluxing_field is
        not specified, it is assumed to be 1.0 everywhere, and the raw flux
        with no local-weighting is returned.

        Additionally, the returned flux is defined as flux *into* the surface,
        not flux *out of* the surface.

        Parameters
        ----------
        field : string
            Any field that can be obtained in a data object.  This is the field
            which will be isocontoured and used as the "local_value" in the
            flux equation.
        value : float
            The value at which the isocontour should be calculated.
        field_x : string
            The x-component field
        field_y : string
            The y-component field
        field_z : string
            The z-component field
        fluxing_field : string, optional
            The field whose passage over the surface is of interest.  If not
            specified, assumed to be 1.0 everywhere.

        Returns
        -------
        flux : float
            The summed flux.  Note that it is not currently scaled; this is
            simply the code-unit area times the fields.

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        calculate the metal flux over it.

        >>> dd = ds.all_data()
        >>> rho = dd.quantities["WeightedAverageQuantity"](
        ...     "Density", weight="CellMassMsun")
        >>> flux = dd.calculate_isocontour_flux("Density", rho,
        ...     "velocity_x", "velocity_y", "velocity_z", "Metal_Density")
        """
        flux = 0.0
        for block, mask in self.blocks:
            flux += self._calculate_flux_in_grid(block, mask, field, value, field_x,
                                                 field_y, field_z, fluxing_field)
        flux = self.comm.mpi_allreduce(flux, op="sum")
        return flux

    def _calculate_flux_in_grid(self, grid, mask, field, value,
                    field_x, field_y, field_z, fluxing_field = None):
        
        vc_fields = [field, field_x, field_y, field_z]
        if fluxing_field is not None:
            vc_fields.append(fluxing_field)

        vc_data = grid.get_vertex_centered_data(vc_fields)

        if fluxing_field is None:
            ff = np.ones_like(vc_data[field], dtype="float64")
        else:
            ff = vc_data[fluxing_field]

        return march_cubes_grid_flux(value, vc_data[field], vc_data[field_x],
            vc_data[field_y], vc_data[field_z], ff, mask, grid.LeftEdge,
            grid.dds)

    def extract_connected_sets(self, field, num_levels, min_val, max_val,
                               log_space=True, cumulative=True):
        """
        This function will create a set of contour objects, defined
        by having connected cell structures, which can then be
        studied and used to 'paint' their source grids, thus enabling
        them to be plotted.

        Note that this function *can* return a connected set object that has no
        member values.
        """
        if log_space:
            cons = np.logspace(np.log10(min_val),np.log10(max_val),
                               num_levels+1)
        else:
            cons = np.linspace(min_val, max_val, num_levels+1)
        contours = {}
        for level in range(num_levels):
            contours[level] = {}
            if cumulative:
                mv = max_val
            else:
                mv = cons[level+1]
            from yt.analysis_modules.level_sets.api import identify_contours
            from yt.analysis_modules.level_sets.clump_handling import \
                add_contour_field
            nj, cids = identify_contours(self, field, cons[level], mv)
            unique_contours = set([])
            for sl_list in cids.values():
                for sl, ff in sl_list:
                    unique_contours.update(np.unique(ff))
            contour_key = uuid.uuid4().hex
            # In case we're a cut region already...
            base_object = getattr(self, 'base_object', self)
            add_contour_field(base_object.ds, contour_key)
            for cid in sorted(unique_contours):
                if cid == -1: continue
                contours[level][cid] = base_object.cut_region(
                    ["obj['contours_%s'] == %s" % (contour_key, cid)],
                    {'contour_slices_%s' % contour_key: cids})
        return cons, contours

    def paint_grids(self, field, value, default_value=None):
        """
        This function paints every cell in our dataset with a given *value*.
        If default_value is given, the other values for the given in every grid
        are discarded and replaced with *default_value*.  Otherwise, the field is
        mandated to 'know how to exist' in the grid.

        Note that this only paints the cells *in the dataset*, so cells in grids
        with child cells are left untouched.
        """
        for grid in self._grids:
            if default_value is not None:
                grid[field] = np.ones(grid.ActiveDimensions)*default_value
            grid[field][self._get_point_indices(grid)] = value

    _particle_handler = None

    @property
    def particles(self):
        if self._particle_handler is None:
            self._particle_handler = \
                particle_handler_registry[self._type_name](self.ds, self)
        return self._particle_handler


    def volume(self):
        """
        Return the volume of the data container.
        This is found by adding up the volume of the cells with centers
        in the container, rather than using the geometric shape of
        the container, so this may vary very slightly
        from what might be expected from the geometric volume.
        """
        return self.quantities.total_quantity(("index", "cell_volume"))

# Many of these items are set up specifically to ensure that
# we are not breaking old pickle files.  This means we must only call the
# _reconstruct_object and that we cannot mandate any additional arguments to
# the reconstruction function.
#
# In the future, this would be better off being set up to more directly
# reference objects or retain state, perhaps with a context manager.
#
# One final detail: time series or multiple datasets in a single pickle
# seems problematic.

class ReconstructedObject(tuple):
    pass

def _check_nested_args(arg, ref_ds):
    if not isinstance(arg, (tuple, list, ReconstructedObject)):
        return arg
    elif isinstance(arg, ReconstructedObject) and ref_ds == arg[0]:
        return arg[1]
    narg = [_check_nested_args(a, ref_ds) for a in arg]
    return narg

def _get_ds_by_hash(hash):
    from yt.data_objects.static_output import _cached_datasets
    for ds in _cached_datasets.values():
        if ds._hash() == hash: return ds
    return None

def _reconstruct_object(*args, **kwargs):
    dsid = args[0]
    dtype = args[1]
    ds = _get_ds_by_hash(dsid)
    if not ds:
        datasets = ParameterFileStore()
        ds = datasets.get_ds_hash(dsid)
    field_parameters = args[-1]
    # will be much nicer when we can do dsid, *a, fp = args
    args = args[2:-1]
    new_args = [_check_nested_args(a, ds) for a in args]
    cls = getattr(ds, dtype)
    obj = cls(*new_args)
    obj.field_parameters.update(field_parameters)
    return ReconstructedObject((ds, obj))

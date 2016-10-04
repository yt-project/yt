"""
Dataset and related data structures.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import functools
import itertools
import numpy as np
import os
import time
import weakref

from collections import defaultdict
from yt.extern.six import add_metaclass, string_types

from yt.config import ytcfg
from yt.fields.derived_field import \
    DerivedField
from yt.funcs import \
    mylog, \
    set_intersection, \
    ensure_list
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.exceptions import \
    YTObjectNotImplemented, \
    YTFieldNotFound, \
    YTGeometryNotSupported
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.utilities.parameter_file_storage import \
    ParameterFileStore, \
    NoParameterShelf, \
    output_type_registry
from yt.units.dimensions import current_mks
from yt.units.unit_object import Unit, unit_system_registry
from yt.units.unit_registry import UnitRegistry
from yt.fields.derived_field import \
    ValidateSpatial
from yt.fields.fluid_fields import \
    setup_gradient_fields
from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field
from yt.data_objects.particle_filters import \
    filter_registry
from yt.data_objects.particle_unions import \
    ParticleUnion
from yt.data_objects.data_containers import \
    data_object_registry
from yt.utilities.minimal_representation import \
    MinimalDataset
from yt.units.yt_array import \
    YTArray, \
    YTQuantity
from yt.units.unit_systems import create_code_unit_system
from yt.data_objects.region_expression import \
    RegionExpression
from yt.geometry.coordinates.api import \
    CoordinateHandler, \
    CartesianCoordinateHandler, \
    PolarCoordinateHandler, \
    CylindricalCoordinateHandler, \
    SphericalCoordinateHandler, \
    GeographicCoordinateHandler, \
    SpectralCubeCoordinateHandler, \
    InternalGeographicCoordinateHandler

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is contant up
# to here, and then have it instantiate EnzoDatasets as appropriate.

_cached_datasets = weakref.WeakValueDictionary()
_ds_store = ParameterFileStore()

def _unsupported_object(ds, obj_name):
    def _raise_unsupp(*args, **kwargs):
        raise YTObjectNotImplemented(ds, obj_name)
    return _raise_unsupp

class RegisteredDataset(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        output_type_registry[name] = cls
        mylog.debug("Registering: %s as %s", name, cls)

class FieldTypeContainer(object):
    def __init__(self, ds):
        self.ds = weakref.proxy(ds)

    def __getattr__(self, attr):
        ds = self.__getattribute__('ds')
        fnc = FieldNameContainer(ds, attr)
        if len(dir(fnc)) == 0:
            return self.__getattribute__(attr)
        return fnc

    def __dir__(self):
        return list(set(t for t, n in self.ds.field_info))

class FieldNameContainer(object):
    def __init__(self, ds, field_type):
        self.ds = ds
        self.field_type = field_type

    def __getattr__(self, attr):
        ft = self.__getattribute__("field_type")
        ds = self.__getattribute__("ds")
        if (ft, attr) not in ds.field_info:
            return self.__getattribute__(attr)
        return ds.field_info[ft, attr]

    def __dir__(self):
        return [n for t, n in self.ds.field_info
                if t == self.field_type]

class IndexProxy(object):
    # This is a simple proxy for Index objects.  It enables backwards
    # compatibility so that operations like .h.sphere, .h.print_stats and
    # .h.grid_left_edge will correctly pass through to the various dataset or
    # index objects.
    def __init__(self, ds):
        self.ds = weakref.proxy(ds)
        ds.index

    def __getattr__(self, name):
        # Check the ds first
        if hasattr(self.ds, name):
            return getattr(self.ds, name)
        # Now for a subset of the available items, check the ds.index.
        elif name in self.ds.index._index_properties:
            return getattr(self.ds.index, name)
        raise AttributeError

class MutableAttribute(object):
    """A descriptor for mutable data"""
    def __init__(self):
        self.data = weakref.WeakKeyDictionary()

    def __get__(self, instance, owner):
        ret = self.data.get(instance, None)
        try:
            ret = ret.copy()
        except AttributeError:
            pass
        return ret

    def __set__(self, instance, value):
        self.data[instance] = value

def requires_index(attr_name):
    @property
    def ireq(self):
        self.index
        # By now it should have been set
        attr = self.__dict__[attr_name]
        return attr

    @ireq.setter
    def ireq(self, value):
        self.__dict__[attr_name] = value

    return ireq

@add_metaclass(RegisteredDataset)
class Dataset(object):

    default_fluid_type = "gas"
    default_field = ("gas", "density")
    fluid_types = ("gas", "deposit", "index")
    particle_types = ("io",) # By default we have an 'all'
    particle_types_raw = ("io",)
    geometry = "cartesian"
    coordinates = None
    max_level = 99
    storage_filename = None
    particle_unions = None
    known_filters = None
    _index_class = None
    field_units = None
    derived_field_list = requires_index("derived_field_list")
    fields = requires_index("fields")
    _instantiated = False
    _particle_type_counts = None

    def __new__(cls, filename=None, *args, **kwargs):
        if not isinstance(filename, string_types):
            obj = object.__new__(cls)
            # The Stream frontend uses a StreamHandler object to pass metadata
            # to __init__.
            is_stream = (hasattr(filename, 'get_fields') and
                         hasattr(filename, 'get_particle_type'))
            if not is_stream:
                obj.__init__(filename, *args, **kwargs)
            return obj
        apath = os.path.abspath(filename)
        if ytcfg.getboolean("yt","skip_dataset_cache"):
            obj = object.__new__(cls)
        elif apath not in _cached_datasets:
            obj = object.__new__(cls)
            if obj._skip_cache is False:
                _cached_datasets[apath] = obj
        else:
            obj = _cached_datasets[apath]
        return obj

    def __init__(self, filename, dataset_type=None, file_style=None, 
                 units_override=None, unit_system="cgs"):
        """
        Base class for generating new output types.  Principally consists of
        a *filename* and a *dataset_type* which will be passed on to children.
        """
        # We return early and do NOT initialize a second time if this file has
        # already been initialized.
        if self._instantiated: return
        self.dataset_type = dataset_type
        self.file_style = file_style
        self.conversion_factors = {}
        self.parameters = {}
        self.region_expression = self.r = RegionExpression(self)
        self.known_filters = self.known_filters or {}
        self.particle_unions = self.particle_unions or {}
        self.field_units = self.field_units or {}
        if units_override is None:
            units_override = {}
        self.units_override = units_override

        # path stuff
        self.parameter_filename = str(filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.expanduser(os.path.dirname(filename))
        self.fullpath = os.path.abspath(self.directory)
        self.backup_filename = self.parameter_filename + '_backup.gdf'
        self.read_from_backup = False
        if os.path.exists(self.backup_filename):
            self.read_from_backup = True
        if len(self.directory) == 0:
            self.directory = "."

        # to get the timing right, do this before the heavy lifting
        self._instantiated = time.time()

        self.min_level = 0
        self.no_cgs_equiv_length = False

        self._create_unit_registry()

        create_code_unit_system(self)
        if unit_system == "code":
            unit_system = str(self)
        else:
            unit_system = str(unit_system).lower()

        self.unit_system = unit_system_registry[unit_system]

        self._parse_parameter_file()
        self.set_units()
        self._setup_coordinate_handler()

        # Because we need an instantiated class to check the ds's existence in
        # the cache, we move that check to here from __new__.  This avoids
        # double-instantiation.
        try:
            _ds_store.check_ds(self)
        except NoParameterShelf:
            pass
        self.print_key_parameters()

        self._set_derived_attrs()
        self._setup_classes()

    def _set_derived_attrs(self):
        if self.domain_left_edge is None or self.domain_right_edge is None:
            self.domain_center = np.zeros(3)
            self.domain_width = np.zeros(3)
        else:
            self.domain_center = 0.5 * (self.domain_right_edge + self.domain_left_edge)
            self.domain_width = self.domain_right_edge - self.domain_left_edge
        if not isinstance(self.current_time, YTQuantity):
            self.current_time = self.quan(self.current_time, "code_time")
        for attr in ("center", "width", "left_edge", "right_edge"):
            n = "domain_%s" % attr
            v = getattr(self, n)
            v = self.arr(v, "code_length")
            setattr(self, n, v)

    def __reduce__(self):
        args = (self._hash(),)
        return (_reconstruct_ds, args)

    def __repr__(self):
        return self.basename

    def _hash(self):
        s = "%s;%s;%s" % (self.basename,
            self.current_time, self.unique_identifier)
        try:
            import hashlib
            return hashlib.md5(s.encode('utf-8')).hexdigest()
        except ImportError:
            return s.replace(";", "*")

    domain_left_edge = MutableAttribute()
    domain_right_edge = MutableAttribute()
    domain_width = MutableAttribute()
    domain_dimensions = MutableAttribute()
    domain_center = MutableAttribute()

    @property
    def _mrep(self):
        return MinimalDataset(self)

    @property
    def _skip_cache(self):
        return False

    def hub_upload(self):
        self._mrep.upload()

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return False

    @classmethod
    def _guess_candidates(cls, base, directories, files):
        """
        This is a class method that accepts a directory (base), a list of files
        in that directory, and a list of subdirectories.  It should return a
        list of filenames (defined relative to the supplied directory) and a
        boolean as to whether or not further directories should be recursed.
        
        This function doesn't need to catch all possibilities, nor does it need
        to filter possibilities.
        """
        return [], True

    def close(self):
        pass

    def __getitem__(self, key):
        """ Returns units, parameters, or conversion_factors in that order. """
        return self.parameters[key]

    def __iter__(self):
      for i in self.parameters: yield i

    def get_smallest_appropriate_unit(self, v, quantity='distance',
                                      return_quantity=False):
        """
        Returns the largest whole unit smaller than the YTQuantity passed to
        it as a string.

        The quantity keyword can be equal to `distance` or `time`.  In the
        case of distance, the units are: 'Mpc', 'kpc', 'pc', 'au', 'rsun',
        'km', etc.  For time, the units are: 'Myr', 'kyr', 'yr', 'day', 'hr',
        's', 'ms', etc.

        If return_quantity is set to True, it finds the largest YTQuantity
        object with a whole unit and a power of ten as the coefficient, and it
        returns this YTQuantity.
        """
        good_u = None
        if quantity == 'distance':
            unit_list =['Ppc', 'Tpc', 'Gpc', 'Mpc', 'kpc', 'pc', 'au', 'rsun',
                        'km', 'cm', 'um', 'nm', 'pm']
        elif quantity == 'time':
            unit_list =['Yyr', 'Zyr', 'Eyr', 'Pyr', 'Tyr', 'Gyr', 'Myr', 'kyr',
                        'yr', 'day', 'hr', 's', 'ms', 'us', 'ns', 'ps', 'fs']
        else:
            raise SyntaxError("Specified quantity must be equal to 'distance'"\
                              "or 'time'.")
        for unit in unit_list:
            uq = self.quan(1.0, unit)
            if uq <= v:
                good_u = unit
                break
        if good_u is None and quantity == 'distance': good_u = 'cm'
        if good_u is None and quantity == 'time': good_u = 's'
        if return_quantity:
            unit_index = unit_list.index(good_u)
            # This avoids indexing errors
            if unit_index == 0: return self.quan(1, unit_list[0])
            # Number of orders of magnitude between unit and next one up
            OOMs = np.ceil(np.log10(self.quan(1, unit_list[unit_index-1]) /
                                    self.quan(1, unit_list[unit_index])))
            # Backwards order of coefficients (e.g. [100, 10, 1])
            coeffs = 10**np.arange(OOMs)[::-1]
            for j in coeffs:
                uq = self.quan(j, good_u)
                if uq <= v:
                    return uq
        else:
            return good_u

    def has_key(self, key):
        """
        Checks units, parameters, and conversion factors. Returns a boolean.

        """
        return key in self.parameters

    _instantiated_index = None
    @property
    def index(self):
        if self._instantiated_index is None:
            if self._index_class is None:
                raise RuntimeError("You should not instantiate Dataset.")
            self._instantiated_index = self._index_class(
                self, dataset_type=self.dataset_type)
            # Now we do things that we need an instantiated index for
            # ...first off, we create our field_info now.
            oldsettings = np.geterr()
            np.seterr(all='ignore')
            self.create_field_info()
            np.seterr(**oldsettings)
        return self._instantiated_index

    _index_proxy = None
    @property
    def h(self):
        if self._index_proxy is None:
            self._index_proxy = IndexProxy(self)
        return self._index_proxy
    hierarchy = h

    @parallel_root_only
    def print_key_parameters(self):
        for a in ["current_time", "domain_dimensions", "domain_left_edge",
                  "domain_right_edge", "cosmological_simulation"]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)
        if hasattr(self, "cosmological_simulation") and \
           getattr(self, "cosmological_simulation"):
            for a in ["current_redshift", "omega_lambda", "omega_matter",
                      "hubble_constant"]:
                if not hasattr(self, a):
                    mylog.error("Missing %s in parameter file definition!", a)
                    continue
                v = getattr(self, a)
                mylog.info("Parameters: %-25s = %s", a, v)

    @parallel_root_only
    def print_stats(self):
        self.index.print_stats()

    @property
    def field_list(self):
        return self.index.field_list

    def create_field_info(self):
        self.field_dependencies = {}
        self.derived_field_list = []
        self.filtered_particle_types = []
        self.field_info = self._field_info_class(self, self.field_list)
        self.coordinates.setup_fields(self.field_info)
        self.field_info.setup_fluid_fields()
        for ptype in self.particle_types:
            self.field_info.setup_particle_fields(ptype)
        self.field_info.setup_fluid_index_fields()
        if "all" not in self.particle_types:
            mylog.debug("Creating Particle Union 'all'")
            pu = ParticleUnion("all", list(self.particle_types_raw))
            nfields = self.add_particle_union(pu)
            if nfields == 0:
                mylog.debug("zero common fields: skipping particle union 'all'")
        self.field_info.setup_extra_union_fields()
        mylog.debug("Loading field plugins.")
        self.field_info.load_all_plugins()
        deps, unloaded = self.field_info.check_derived_fields()
        self.field_dependencies.update(deps)
        self.fields = FieldTypeContainer(self)
        self.index.field_list = sorted(self.field_list)

    def setup_deprecated_fields(self):
        from yt.fields.field_aliases import _field_name_aliases
        added = []
        for old_name, new_name in _field_name_aliases:
            try:
                fi = self._get_field_info(new_name)
            except YTFieldNotFound:
                continue
            self.field_info.alias(("gas", old_name), fi.name)
            added.append(("gas", old_name))
        self.field_info.find_dependencies(added)

    def _setup_coordinate_handler(self):
        kwargs = {}
        if isinstance(self.geometry, tuple):
            self.geometry, ordering = self.geometry
            kwargs['ordering'] = ordering
        if isinstance(self.geometry, CoordinateHandler):
            # I kind of dislike this.  The geometry field should always be a
            # string, but the way we're set up with subclassing, we can't
            # mandate that quite the way I'd like.
            self.coordinates = self.geometry
            return
        elif callable(self.geometry):
            cls = self.geometry
        elif self.geometry == "cartesian":
            cls = CartesianCoordinateHandler
        elif self.geometry == "cylindrical":
            cls = CylindricalCoordinateHandler
        elif self.geometry == "polar":
            cls = PolarCoordinateHandler
        elif self.geometry == "spherical":
            cls = SphericalCoordinateHandler
        elif self.geometry == "geographic":
            cls = GeographicCoordinateHandler
        elif self.geometry == "internal_geographic":
            cls = InternalGeographicCoordinateHandler
        elif self.geometry == "spectral_cube":
            cls = SpectralCubeCoordinateHandler
        else:
            raise YTGeometryNotSupported(self.geometry)
        self.coordinates = cls(self, **kwargs)

    def add_particle_union(self, union):
        # No string lookups here, we need an actual union.
        f = self.particle_fields_by_type

        # find fields common to all particle types in the union
        fields = set_intersection(
            [f[s] for s in union if s in self.particle_types_raw]
        )

        if len(fields) == 0:
            # don't create this union if no fields are common to all
            # particle types
            return len(fields)

        for field in fields:
            units = set([])
            for s in union:
                # First we check our existing fields for units
                funits = self._get_field_info(s, field).units
                # Then we override with field_units settings.
                funits = self.field_units.get((s, field), funits)
                units.add(funits)
            if len(units) == 1:
                self.field_units[union.name, field] = list(units)[0]
        self.particle_types += (union.name,)
        self.particle_unions[union.name] = union
        fields = [ (union.name, field) for field in fields]
        new_fields = [_ for _ in fields if _ not in self.field_list]
        self.field_list.extend(new_fields)
        new_field_info_fields = [
            _ for _ in fields if _ not in self.field_info.field_list]
        self.field_info.field_list.extend(new_field_info_fields)
        self.index.field_list = sorted(self.field_list)
        # Give ourselves a chance to add them here, first, then...
        # ...if we can't find them, we set them up as defaults.
        new_fields = self._setup_particle_types([union.name])
        self.field_info.find_dependencies(new_fields)
        return len(new_fields)

    def add_particle_filter(self, filter):
        # This requires an index
        self.index
        # This is a dummy, which we set up to enable passthrough of "all"
        # concatenation fields.
        n = getattr(filter, "name", filter)
        self.known_filters[n] = None
        if isinstance(filter, string_types):
            used = False
            for f in filter_registry[filter]:
                used = self._setup_filtered_type(f)
                if used:
                    filter = f
                    break
        else:
            used = self._setup_filtered_type(filter)
        if not used:
            self.known_filters.pop(n, None)
            return False
        self.known_filters[filter.name] = filter
        return True

    def _setup_filtered_type(self, filter):
        if not filter.available(self.derived_field_list):
            return False
        fi = self.field_info
        fd = self.field_dependencies
        available = False
        for fn in self.derived_field_list:
            if fn[0] == filter.filtered_type:
                # Now we can add this
                available = True
                self.derived_field_list.append(
                    (filter.name, fn[1]))
                fi[filter.name, fn[1]] = filter.wrap_func(fn, fi[fn])
                # Now we append the dependencies
                fd[filter.name, fn[1]] = fd[fn]
        if available:
            self.particle_types += (filter.name,)
            self.filtered_particle_types.append(filter.name)
            new_fields = self._setup_particle_types([filter.name])
            deps, _ = self.field_info.check_derived_fields(new_fields)
            self.field_dependencies.update(deps)
        return available

    def _setup_particle_types(self, ptypes = None):
        df = []
        if ptypes is None: ptypes = self.ds.particle_types_raw
        for ptype in set(ptypes):
            df += self._setup_particle_type(ptype)
        return df

    _last_freq = (None, None)
    _last_finfo = None
    def _get_field_info(self, ftype, fname = None):
        self.index
        if fname is None:
            if isinstance(ftype, DerivedField):
                ftype, fname = ftype.name
            else:
                ftype, fname = "unknown", ftype
        guessing_type = False
        if ftype == "unknown":
            guessing_type = True
            ftype = self._last_freq[0] or ftype
        field = (ftype, fname)
        if field == self._last_freq:
            if field not in self.field_info.field_aliases.values():
                return self._last_finfo
        if field in self.field_info:
            self._last_freq = field
            self._last_finfo = self.field_info[(ftype, fname)]
            return self._last_finfo
        if fname in self.field_info:
            # Sometimes, if guessing_type == True, this will be switched for
            # the type of field it is.  So we look at the field type and
            # determine if we need to change the type.
            fi = self._last_finfo = self.field_info[fname]
            if fi.particle_type and self._last_freq[0] \
                not in self.particle_types:
                    field = "all", field[1]
            elif not fi.particle_type and self._last_freq[0] \
                not in self.fluid_types:
                    field = self.default_fluid_type, field[1]
            self._last_freq = field
            return self._last_finfo
        # We also should check "all" for particles, which can show up if you're
        # mixing deposition/gas fields with particle fields.
        if guessing_type:
            to_guess = ["all", self.default_fluid_type] \
                     + list(self.fluid_types) \
                     + list(self.particle_types)
            for ftype in to_guess:
                if (ftype, fname) in self.field_info:
                    self._last_freq = (ftype, fname)
                    self._last_finfo = self.field_info[(ftype, fname)]
                    return self._last_finfo
        raise YTFieldNotFound((ftype, fname), self)

    def _setup_classes(self):
        # Called by subclass
        self.object_types = []
        self.objects = []
        self.plots = []
        for name, cls in sorted(data_object_registry.items()):
            if name in self._index_class._unsupported_objects:
                setattr(self, name,
                    _unsupported_object(self, name))
                continue
            cname = cls.__name__
            if cname.endswith("Base"): cname = cname[:-4]
            self._add_object_class(name, cls)
        if self.refine_by != 2 and hasattr(self, 'proj') and \
            hasattr(self, 'overlap_proj'):
            mylog.warning("Refine by something other than two: reverting to"
                        + " overlap_proj")
            self.proj = self.overlap_proj
        if self.dimensionality < 3 and hasattr(self, 'proj') and \
            hasattr(self, 'overlap_proj'):
            mylog.warning("Dimensionality less than 3: reverting to"
                        + " overlap_proj")
            self.proj = self.overlap_proj
        self.object_types.sort()

    def _add_object_class(self, name, base):
        self.object_types.append(name)
        obj = functools.partial(base, ds=weakref.proxy(self))
        obj.__doc__ = base.__doc__
        setattr(self, name, obj)

    def find_max(self, field):
        """
        Returns (value, location) of the maximum of a given field.
        """
        mylog.debug("Searching for maximum value of %s", field)
        source = self.all_data()
        max_val, mx, my, mz = \
            source.quantities.max_location(field)
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f",
              max_val, mx, my, mz)
        return max_val, self.arr([mx, my, mz], 'code_length', dtype="float64")

    def find_min(self, field):
        """
        Returns (value, location) for the minimum of a given field.
        """
        mylog.debug("Searching for minimum value of %s", field)
        source = self.all_data()
        min_val, mx, my, mz = \
            source.quantities.min_location(field)
        mylog.info("Min Value is %0.5e at %0.16f %0.16f %0.16f",
              min_val, mx, my, mz)
        return min_val, self.arr([mx, my, mz], 'code_length', dtype="float64")

    def find_field_values_at_point(self, fields, coords):
        """
        Returns the values [field1, field2,...] of the fields at the given
        coordinates. Returns a list of field values in the same order as
        the input *fields*.
        """
        point = self.point(coords)
        ret = []
        field_list = ensure_list(fields)
        for field in field_list:
            ret.append(point[field])
        if len(field_list) == 1:
            return ret[0]
        else:
            return ret

    def find_field_values_at_points(self, fields, coords):
        """
        Returns the values [field1, field2,...] of the fields at the given
        [(x1, y1, z2), (x2, y2, z2),...] points.  Returns a list of field
        values in the same order as the input *fields*.

        """
        # If an optimized version exists on the Index object we'll use that
        try:
            return self.index._find_field_values_at_points(fields, coords)
        except AttributeError:
            pass

        fields = ensure_list(fields)
        out = []

        # This may be slow because it creates a data object for each point
        for field_index, field in enumerate(fields):
            funit = self._get_field_info[field].units
            out.append(self.arr(np.empty((len(coords),)), funit))
            for coord_index, coord in enumerate(coords):
                out[field_index][coord_index] = self.point(coord)[fields]
        if len(fields) == 1:
            return out[0]
        else:
            return out

    # Now all the object related stuff
    def all_data(self, find_max=False, **kwargs):
        """
        all_data is a wrapper to the Region object for creating a region
        which covers the entire simulation domain.
        """
        if find_max: c = self.find_max("density")[1]
        else: c = (self.domain_right_edge + self.domain_left_edge)/2.0
        return self.region(c,
            self.domain_left_edge, self.domain_right_edge, **kwargs)

    def box(self, left_edge, right_edge, **kwargs):
        """
        box is a wrapper to the Region object for creating a region
        without having to specify a *center* value.  It assumes the center
        is the midpoint between the left_edge and right_edge.
        """
        # we handle units in the region data object
        # but need to check if left_edge or right_edge is a
        # list or other non-array iterable before calculating
        # the center
        if not isinstance(left_edge, np.ndarray):
            left_edge = np.array(left_edge)
        if not isinstance(right_edge, np.ndarray):
            right_edge = np.array(right_edge)
        c = (left_edge + right_edge)/2.0
        return self.region(c, left_edge, right_edge, **kwargs)

    def _setup_particle_type(self, ptype):
        orig = set(self.field_info.items())
        self.field_info.setup_particle_fields(ptype)
        return [n for n, v in set(self.field_info.items()).difference(orig)]

    @property
    def particle_fields_by_type(self):
        fields = defaultdict(list)
        for field in self.field_list:
            if field[0] in self.particle_types_raw:
                fields[field[0]].append(field[1])
        return fields

    @property
    def particles_exist(self):
        for pt, f in itertools.product(self.particle_types_raw, self.field_list):
            if pt == f[0]:
                return True
        return False

    @property
    def particle_type_counts(self):
        self.index
        if self.particles_exist is False:
            return {}

        # frontends or index implementation can populate this dict while
        # creating the index if they know particle counts at that time
        if self._particle_type_counts is not None:
            return self._particle_type_counts

        self._particle_type_counts = self.index._get_particle_type_counts()
        return self._particle_type_counts

    @property
    def ires_factor(self):
        o2 = np.log2(self.refine_by)
        if o2 != int(o2):
            raise RuntimeError
        return int(o2)

    def relative_refinement(self, l0, l1):
        return self.refine_by**(l1-l0)

    def _create_unit_registry(self):
        self.unit_registry = UnitRegistry()
        import yt.units.dimensions as dimensions
        self.unit_registry.add("code_length", 1.0, dimensions.length)
        self.unit_registry.add("code_mass", 1.0, dimensions.mass)
        self.unit_registry.add("code_density", 1.0, dimensions.density)
        self.unit_registry.add("code_time", 1.0, dimensions.time)
        self.unit_registry.add("code_magnetic", 1.0, dimensions.magnetic_field)
        self.unit_registry.add("code_temperature", 1.0, dimensions.temperature)
        self.unit_registry.add("code_pressure", 1.0, dimensions.pressure)
        self.unit_registry.add("code_velocity", 1.0, dimensions.velocity)
        self.unit_registry.add("code_metallicity", 1.0,
                               dimensions.dimensionless)
        self.unit_registry.add("a", 1.0, dimensions.dimensionless)

    def set_units(self):
        """
        Creates the unit registry for this dataset.

        """
        from yt.units.dimensions import length
        if getattr(self, "cosmological_simulation", False):
            # this dataset is cosmological, so add cosmological units.
            self.unit_registry.modify("h", self.hubble_constant)
            # Comoving lengths
            for my_unit in ["m", "pc", "AU", "au"]:
                new_unit = "%scm" % my_unit
                self.unit_registry.add(new_unit, self.unit_registry.lut[my_unit][0] /
                                       (1 + self.current_redshift),
                                       length, "\\rm{%s}/(1+z)" % my_unit)
            self.unit_registry.modify('a', 1/(1+self.current_redshift))

        self.set_code_units()

        if getattr(self, "cosmological_simulation", False):
            # this dataset is cosmological, add a cosmology object
            self.cosmology = \
                    Cosmology(hubble_constant=self.hubble_constant,
                              omega_matter=self.omega_matter,
                              omega_lambda=self.omega_lambda)
            self.critical_density = \
                    self.cosmology.critical_density(self.current_redshift)
            self.scale_factor = 1.0 / (1.0 + self.current_redshift)

    def get_unit_from_registry(self, unit_str):
        """
        Creates a unit object matching the string expression, using this
        dataset's unit registry.

        Parameters
        ----------
        unit_str : str
            string that we can parse for a sympy Expr.

        """
        new_unit = Unit(unit_str, registry=self.unit_registry)
        return new_unit

    def set_code_units(self):
        # here we override units, if overrides have been provided.
        self._override_code_units()

        # set attributes like ds.length_unit
        self._set_code_unit_attributes()

        self.unit_registry.modify("code_length", self.length_unit)
        self.unit_registry.modify("code_mass", self.mass_unit)
        self.unit_registry.modify("code_time", self.time_unit)
        if hasattr(self, 'magnetic_unit'):
            # if the magnetic unit is in T, we need to recreate the code unit
            # system as an MKS-like system
            if current_mks in self.magnetic_unit.units.dimensions.free_symbols:

                if self.unit_system == unit_system_registry[str(self)]:
                    unit_system_registry.pop(str(self))
                    create_code_unit_system(self, current_mks_unit='A')
                    self.unit_system = unit_system_registry[str(self)]
                elif str(self.unit_system) == 'mks':
                    pass
                else:
                    self.magnetic_unit = \
                        self.magnetic_unit.to_equivalent('gauss', 'CGS')
            self.unit_registry.modify("code_magnetic", self.magnetic_unit)
        vel_unit = getattr(
            self, "velocity_unit", self.length_unit / self.time_unit)
        pressure_unit = getattr(
            self, "pressure_unit",
            self.mass_unit / (self.length_unit * (self.time_unit)**2))
        temperature_unit = getattr(self, "temperature_unit", 1.0)
        density_unit = getattr(self, "density_unit", self.mass_unit / self.length_unit**3)
        self.unit_registry.modify("code_velocity", vel_unit)
        self.unit_registry.modify("code_temperature", temperature_unit)
        self.unit_registry.modify("code_pressure", pressure_unit)
        self.unit_registry.modify("code_density", density_unit)
        # domain_width does not yet exist
        if (self.domain_left_edge is not None and
            self.domain_right_edge is not None):
            DW = self.arr(self.domain_right_edge - self.domain_left_edge, "code_length")
            self.unit_registry.add("unitary", float(DW.max() * DW.units.base_value),
                                   DW.units.dimensions)

    def _override_code_units(self):
        if len(self.units_override) == 0:
            return
        mylog.warning(
            "Overriding code units. This is an experimental and potentially "
            "dangerous option that may yield inconsistent results, and must be "
            "used very carefully, and only if you know what you want from it.")
        for unit, cgs in [("length", "cm"), ("time", "s"), ("mass", "g"),
                          ("velocity","cm/s"), ("magnetic","gauss"), 
                          ("temperature","K")]:
            val = self.units_override.get("%s_unit" % unit, None)
            if val is not None:
                if isinstance(val, YTQuantity):
                    val = (val.v, str(val.units))
                elif not isinstance(val, tuple):
                    val = (val, cgs)
                mylog.info("Overriding %s_unit: %g %s.", unit, val[0], val[1])
                setattr(self, "%s_unit" % unit, self.quan(val[0], val[1]))

    _arr = None
    @property
    def arr(self):
        """Converts an array into a :class:`yt.units.yt_array.YTArray`

        The returned YTArray will be dimensionless by default, but can be
        cast to arbitrary units using the ``input_units`` keyword argument.

        Parameters
        ----------

        input_array : iterable
            A tuple, list, or array to attach units to
        input_units : String unit specification, unit symbol or astropy object
            The units of the array. Powers must be specified using python syntax
            (cm**3, not cm^3).
        dtype : string or NumPy dtype object
            The dtype of the returned array data

        Examples
        --------

        >>> import yt
        >>> import numpy as np
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> a = ds.arr([1, 2, 3], 'cm')
        >>> b = ds.arr([4, 5, 6], 'm')
        >>> a + b
        YTArray([ 401.,  502.,  603.]) cm
        >>> b + a
        YTArray([ 4.01,  5.02,  6.03]) m

        Arrays returned by this function know about the dataset's unit system

        >>> a = ds.arr(np.ones(5), 'code_length')
        >>> a.in_units('Mpccm/h')
        YTArray([ 1.00010449,  1.00010449,  1.00010449,  1.00010449,
                 1.00010449]) Mpc

        """

        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(YTArray, registry = self.unit_registry)
        return self._arr

    _quan = None
    @property
    def quan(self):
        """Converts an scalar into a :class:`yt.units.yt_array.YTQuantity`

        The returned YTQuantity will be dimensionless by default, but can be
        cast to arbitrary units using the ``input_units`` keyword argument.

        Parameters
        ----------

        input_scalar : an integer or floating point scalar
            The scalar to attach units to
        input_units : String unit specification, unit symbol or astropy object
            The units of the quantity. Powers must be specified using python
            syntax (cm**3, not cm^3).
        dtype : string or NumPy dtype object
            The dtype of the array data.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')

        >>> a = ds.quan(1, 'cm')
        >>> b = ds.quan(2, 'm')
        >>> a + b
        201.0 cm
        >>> b + a
        2.01 m

        Quantities created this way automatically know about the unit system
        of the dataset.

        >>> a = ds.quan(5, 'code_length')
        >>> a.in_cgs()
        1.543e+25 cm

        """

        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(YTQuantity, registry=self.unit_registry)
        return self._quan

    def add_field(self, name, function=None, **kwargs):
        """
        Dataset-specific call to add_field

        Add a new field, along with supplemental metadata, to the list of
        available fields.  This respects a number of arguments, all of which
        are passed on to the constructor for
        :class:`~yt.data_objects.api.DerivedField`.

        Parameters
        ----------

        name : str
           is the name of the field.
        function : callable
           A function handle that defines the field.  Should accept
           arguments (field, data)
        units : str
           A plain text string encoding the unit.  Powers must be in
           python syntax (** instead of ^).
        take_log : bool
           Describes whether the field should be logged
        validators : list
           A list of :class:`FieldValidator` objects
        particle_type : bool
           Is this a particle (1D) field?
        vector_field : bool
           Describes the dimensionality of the field.  Currently unused.
        display_name : str
           A name used in the plots

        """
        self.index
        override = kwargs.get("force_override", False)
        # Handle the case where the field has already been added.
        if not override and name in self.field_info:
            mylog.warning("Field %s already exists. To override use " +
                          "force_override=True.", name)
        self.field_info.add_field(name, function=function, **kwargs)
        self.field_info._show_field_errors.append(name)
        deps, _ = self.field_info.check_derived_fields([name])
        self.field_dependencies.update(deps)

    def add_deposited_particle_field(self, deposit_field, method, kernel_name='cubic',
                                     weight_field='particle_mass'):
        """Add a new deposited particle field

        Creates a new deposited field based on the particle *deposit_field*.

        Parameters
        ----------

        deposit_field : tuple
           The field name tuple of the particle field the deposited field will
           be created from.  This must be a field name tuple so yt can
           appropriately infer the correct particle type.
        method : string
           This is the "method name" which will be looked up in the
           `particle_deposit` namespace as `methodname_deposit`.  Current
           methods include `simple_smooth`, `sum`, `std`, `cic`, `weighted_mean`,
           `mesh_id`, and `nearest`.
        kernel_name : string, default 'cubic'
           This is the name of the smoothing kernel to use. It is only used for
           the `simple_smooth` method and is otherwise ignored. Current
           supported kernel names include `cubic`, `quartic`, `quintic`,
           `wendland2`, `wendland4`, and `wendland6`.
        weight_field : string, default 'particle_mass'
           Weighting field name for deposition method `weighted_mean`.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        self.index
        if isinstance(deposit_field, tuple):
            ptype, deposit_field = deposit_field[0], deposit_field[1]
        else:
            raise RuntimeError

        units = self.field_info[ptype, deposit_field].units
        take_log = self.field_info[ptype, deposit_field].take_log
        name_map = {"sum": "sum", "std":"std", "cic": "cic", "weighted_mean": "avg",
                    "nearest": "nn", "simple_smooth": "ss", "count": "count"}
        field_name = "%s_" + name_map[method] + "_%s"
        field_name = field_name % (ptype, deposit_field.replace('particle_', ''))

        if method == "count":
            field_name = "%s_count" % ptype
            if ("deposit", field_name) in self.field_info:
                mylog.warning("The deposited field %s already exists" % field_name)
                return ("deposit", field_name)
            else:
                units = "dimensionless"
                take_log = False

        def _deposit_field(field, data):
            """
            Create a grid field for particle quantities using given method.
            """
            pos = data[ptype, "particle_position"]
            if method == 'weighted_mean':
                d = data.ds.arr(data.deposit(pos, [data[ptype, deposit_field],
                                                   data[ptype, weight_field]],
                                             method=method, kernel_name=kernel_name),
                                             input_units=units)
                d[np.isnan(d)] = 0.0
            else:
                d = data.ds.arr(data.deposit(pos, [data[ptype, deposit_field]],
                                             method=method, kernel_name=kernel_name),
                                             input_units=units)
            return d

        self.add_field(
            ("deposit", field_name),
            function=_deposit_field,
            units=units,
            take_log=take_log,
            validators=[ValidateSpatial()])
        return ("deposit", field_name)

    def add_smoothed_particle_field(self, smooth_field, method="volume_weighted",
                                    nneighbors=64, kernel_name="cubic"):
        """Add a new smoothed particle field

        Creates a new smoothed field based on the particle *smooth_field*.

        Parameters
        ----------

        smooth_field : tuple
           The field name tuple of the particle field the smoothed field will
           be created from.  This must be a field name tuple so yt can
           appropriately infer the correct particle type.
        method : string, default 'volume_weighted'
           The particle smoothing method to use. Can only be 'volume_weighted'
           for now.
        nneighbors : int, default 64
            The number of neighbors to examine during the process.
        kernel_name : string, default 'cubic'
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        self.index
        if isinstance(smooth_field, tuple):
            ptype, smooth_field = smooth_field[0], smooth_field[1]
        else:
            raise RuntimeError("smooth_field must be a tuple, received %s" %
                               smooth_field)
        if method != "volume_weighted":
            raise NotImplementedError("method must be 'volume_weighted'")

        coord_name = "particle_position"
        mass_name = "particle_mass"
        smoothing_length_name = "smoothing_length"
        if (ptype, smoothing_length_name) not in self.derived_field_list:
            raise ValueError("%s not in derived_field_list" %
                             ((ptype, smoothing_length_name),))
        density_name = "density"
        registry = self.field_info

        return add_volume_weighted_smoothed_field(ptype, coord_name, mass_name,
                   smoothing_length_name, density_name, smooth_field, registry,
                   nneighbors=nneighbors, kernel_name=kernel_name)[0]

    def add_gradient_fields(self, input_field):
        """Add gradient fields.

        Creates four new grid-based fields that represent the components of
        the gradient of an existing field, plus an extra field for the magnitude
        of the gradient. Currently only supported in Cartesian geometries. The
        gradient is computed using second-order centered differences.

        Parameters
        ----------
        input_field : tuple
           The field name tuple of the particle field the deposited field will
           be created from.  This must be a field name tuple so yt can
           appropriately infer the correct field type.

        Returns
        -------
        A list of field name tuples for the newly created fields.

        Examples
        --------
        >>> grad_fields = ds.add_gradient_fields(("gas","temperature"))
        >>> print(grad_fields)
        [('gas', 'temperature_gradient_x'),
         ('gas', 'temperature_gradient_y'),
         ('gas', 'temperature_gradient_z'),
         ('gas', 'temperature_gradient_magnitude')]
        """
        self.index
        if isinstance(input_field, tuple):
            ftype, input_field = input_field[0], input_field[1]
        else:
            raise RuntimeError
        units = self.field_info[ftype, input_field].units
        setup_gradient_fields(self.field_info, (ftype, input_field), units)
        # Now we make a list of the fields that were just made, to check them
        # and to return them
        grad_fields = [(ftype,input_field+"_gradient_%s" % suffix)
                       for suffix in "xyz"]
        grad_fields.append((ftype,input_field+"_gradient_magnitude"))
        deps, _ = self.field_info.check_derived_fields(grad_fields)
        self.field_dependencies.update(deps)
        return grad_fields

def _reconstruct_ds(*args, **kwargs):
    datasets = ParameterFileStore()
    ds = datasets.get_ds_hash(*args)
    return ds

class ParticleFile(object):
    def __init__(self, ds, io, filename, file_id):
        self.ds = ds
        self.io = weakref.proxy(io)
        self.filename = filename
        self.file_id = file_id
        self.total_particles = self.io._count_particles(self)

    def select(self, selector):
        pass

    def count(self, selector):
        pass

    def _calculate_offsets(self, fields):
        pass

    def __lt__(self, other):
        return self.filename < other.filename

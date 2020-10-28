import abc
import functools
import itertools
import os
import pickle
import time
import weakref
from collections import defaultdict
from stat import ST_CTIME

import numpy as np
from unyt.exceptions import UnitConversionError, UnitParseError

from yt.config import ytcfg
from yt.data_objects.particle_filters import filter_registry
from yt.data_objects.particle_unions import ParticleUnion
from yt.data_objects.region_expression import RegionExpression
from yt.fields.derived_field import ValidateSpatial
from yt.fields.field_type_container import FieldTypeContainer
from yt.fields.fluid_fields import setup_gradient_fields
from yt.fields.particle_fields import DEP_MSG_SMOOTH_FIELD
from yt.funcs import (
    ensure_list,
    issue_deprecation_warning,
    iterable,
    mylog,
    set_intersection,
    setdefaultattr,
)
from yt.geometry.coordinates.api import (
    CartesianCoordinateHandler,
    CoordinateHandler,
    CylindricalCoordinateHandler,
    GeographicCoordinateHandler,
    InternalGeographicCoordinateHandler,
    PolarCoordinateHandler,
    SpectralCubeCoordinateHandler,
    SphericalCoordinateHandler,
)
from yt.units import UnitContainer, _wrap_display_ytarray, dimensions
from yt.units.dimensions import current_mks
from yt.units.unit_object import Unit, define_unit
from yt.units.unit_registry import UnitRegistry
from yt.units.unit_systems import create_code_unit_system, unit_system_registry
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.cosmology import Cosmology
from yt.utilities.exceptions import (
    YTFieldNotFound,
    YTGeometryNotSupported,
    YTIllDefinedParticleFilter,
    YTObjectNotImplemented,
)
from yt.utilities.minimal_representation import MinimalDataset
from yt.utilities.object_registries import data_object_registry, output_type_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_root_only
from yt.utilities.parameter_file_storage import NoParameterShelf, ParameterFileStore

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is constant up
# to here, and then have it instantiate EnzoDatasets as appropriate.

_cached_datasets = weakref.WeakValueDictionary()
_ds_store = ParameterFileStore()


def _unsupported_object(ds, obj_name):
    def _raise_unsupp(*args, **kwargs):
        raise YTObjectNotImplemented(ds, obj_name)

    return _raise_unsupp


class IndexProxy:
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


class MutableAttribute:
    """A descriptor for mutable data"""

    def __init__(self, display_array=False):
        self.data = weakref.WeakKeyDictionary()
        self.display_array = display_array

    def __get__(self, instance, owner):
        if not instance:
            return None
        ret = self.data.get(instance, None)
        try:
            ret = ret.copy()
        except AttributeError:
            pass
        if self.display_array:
            try:
                ret._ipython_display_ = functools.partial(_wrap_display_ytarray, ret)
            # This will error out if the items have yet to be turned into
            # YTArrays, in which case we just let it go.
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


class Dataset(abc.ABC):

    default_fluid_type = "gas"
    default_field = ("gas", "density")
    fluid_types = ("gas", "deposit", "index")
    particle_types = ("io",)  # By default we have an 'all'
    particle_types_raw = ("io",)
    geometry = "cartesian"
    coordinates = None
    storage_filename = None
    particle_unions = None
    known_filters = None
    _index_class = None
    field_units = None
    derived_field_list = requires_index("derived_field_list")
    fields = requires_index("fields")
    _instantiated = False
    _unique_identifier = None
    _particle_type_counts = None
    _proj_type = "quad_proj"
    _ionization_label_format = "roman_numeral"

    # these are set in self._parse_parameter_file()
    domain_left_edge = MutableAttribute(True)
    domain_right_edge = MutableAttribute(True)
    domain_dimensions = MutableAttribute(True)
    periodicity = MutableAttribute()

    # these are set in self._set_derived_attrs()
    domain_width = MutableAttribute(True)
    domain_center = MutableAttribute(True)

    def __new__(cls, filename=None, *args, **kwargs):
        if not isinstance(filename, str):
            obj = object.__new__(cls)
            # The Stream frontend uses a StreamHandler object to pass metadata
            # to __init__.
            is_stream = hasattr(filename, "get_fields") and hasattr(
                filename, "get_particle_type"
            )
            if not is_stream:
                obj.__init__(filename, *args, **kwargs)
            return obj
        apath = os.path.abspath(filename)
        cache_key = (apath, pickle.dumps(args), pickle.dumps(kwargs))
        if ytcfg.getboolean("yt", "skip_dataset_cache"):
            obj = object.__new__(cls)
        elif cache_key not in _cached_datasets:
            obj = object.__new__(cls)
            if not obj._skip_cache:
                _cached_datasets[cache_key] = obj
        else:
            obj = _cached_datasets[cache_key]
        return obj

    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(*args, **kwargs)
        output_type_registry[cls.__name__] = cls
        mylog.debug("Registering: %s as %s", cls.__name__, cls)

    def __init__(
        self,
        filename,
        dataset_type=None,
        file_style=None,
        units_override=None,
        unit_system="cgs",
    ):
        """
        Base class for generating new output types.  Principally consists of
        a *filename* and a *dataset_type* which will be passed on to children.
        """
        # We return early and do NOT initialize a second time if this file has
        # already been initialized.
        if self._instantiated:
            return
        self.dataset_type = dataset_type
        self.file_style = file_style
        self.conversion_factors = {}
        self.parameters = {}
        self.region_expression = self.r = RegionExpression(self)
        self.known_filters = self.known_filters or {}
        self.particle_unions = self.particle_unions or {}
        self.field_units = self.field_units or {}
        self.units_override = self.__class__._sanitize_units_override(units_override)

        # path stuff
        self.parameter_filename = str(filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.expanduser(os.path.dirname(filename))
        self.fullpath = os.path.abspath(self.directory)
        self.backup_filename = self.parameter_filename + "_backup.gdf"
        self.read_from_backup = False
        if os.path.exists(self.backup_filename):
            self.read_from_backup = True
        if len(self.directory) == 0:
            self.directory = "."

        # to get the timing right, do this before the heavy lifting
        self._instantiated = time.time()

        self.no_cgs_equiv_length = False

        if unit_system == "code":
            # create a fake MKS unit system which we will override later to
            # avoid chicken/egg issue of the unit registry needing a unit system
            # but code units need a unit registry to define the code units on
            used_unit_system = "mks"
        else:
            used_unit_system = unit_system

        self._create_unit_registry(used_unit_system)

        self._parse_parameter_file()
        self.set_units()
        self.setup_cosmology()
        self._assign_unit_system(unit_system)
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

    @property
    def unique_identifier(self):
        if self._unique_identifier is None:
            self._unique_identifier = int(os.stat(self.parameter_filename)[ST_CTIME])
        return self._unique_identifier

    @unique_identifier.setter
    def unique_identifier(self, value):
        self._unique_identifier = value

    # abstract methods require implementation in subclasses
    @classmethod
    @abc.abstractmethod
    def _is_valid(cls, *args, **kwargs):
        # A heuristic test to determine if the data format can be interpreted
        # with the present frontend
        return False

    @abc.abstractmethod
    def _parse_parameter_file(self):
        # set up various attributes from self.parameter_filename
        # for a full description of what is required here see
        # yt.frontends._skeleton.SkeletonDataset
        pass

    @abc.abstractmethod
    def _set_code_unit_attributes(self):
        # set up code-units to physical units normalization factors
        # for a full description of what is required here see
        # yt.frontends._skeleton.SkeletonDataset
        pass

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
            n = f"domain_{attr}"
            v = getattr(self, n)
            if not isinstance(v, YTArray) and v is not None:
                # Note that we don't add on _ipython_display_ here because
                # everything is stored inside a MutableAttribute.
                v = self.arr(v, "code_length")
                setattr(self, n, v)

    def __reduce__(self):
        args = (self._hash(),)
        return (_reconstruct_ds, args)

    def __repr__(self):
        return self.basename

    def _hash(self):
        s = f"{self.basename};{self.current_time};{self.unique_identifier}"
        try:
            import hashlib

            return hashlib.md5(s.encode("utf-8")).hexdigest()
        except ImportError:
            return s.replace(";", "*")

    _checksum = None

    @property
    def checksum(self):
        """
        Computes md5 sum of a dataset.

        Note: Currently this property is unable to determine a complete set of
        files that are a part of a given dataset. As a first approximation, the
        checksum of :py:attr:`~parameter_file` is calculated. In case
        :py:attr:`~parameter_file` is a directory, checksum of all files inside
        the directory is calculated.
        """
        if self._checksum is None:
            try:
                import hashlib
            except ImportError:
                self._checksum = "nohashlib"
                return self._checksum

            def generate_file_md5(m, filename, blocksize=2 ** 20):
                with open(filename, "rb") as f:
                    while True:
                        buf = f.read(blocksize)
                        if not buf:
                            break
                        m.update(buf)

            m = hashlib.md5()
            if os.path.isdir(self.parameter_filename):
                for root, _, files in os.walk(self.parameter_filename):
                    for fname in files:
                        fname = os.path.join(root, fname)
                        generate_file_md5(m, fname)
            elif os.path.isfile(self.parameter_filename):
                generate_file_md5(m, self.parameter_filename)
            else:
                m = "notafile"

            if hasattr(m, "hexdigest"):
                m = m.hexdigest()
            self._checksum = m
        return self._checksum

    @property
    def _mrep(self):
        return MinimalDataset(self)

    @property
    def _skip_cache(self):
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
        for i in self.parameters:
            yield i

    def get_smallest_appropriate_unit(
        self, v, quantity="distance", return_quantity=False
    ):
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
        if quantity == "distance":
            unit_list = [
                "Ppc",
                "Tpc",
                "Gpc",
                "Mpc",
                "kpc",
                "pc",
                "au",
                "rsun",
                "km",
                "cm",
                "um",
                "nm",
                "pm",
            ]
        elif quantity == "time":
            unit_list = [
                "Yyr",
                "Zyr",
                "Eyr",
                "Pyr",
                "Tyr",
                "Gyr",
                "Myr",
                "kyr",
                "yr",
                "day",
                "hr",
                "s",
                "ms",
                "us",
                "ns",
                "ps",
                "fs",
            ]
        else:
            raise SyntaxError(
                "Specified quantity must be equal to 'distance'" "or 'time'."
            )
        for unit in unit_list:
            uq = self.quan(1.0, unit)
            if uq <= v:
                good_u = unit
                break
        if good_u is None and quantity == "distance":
            good_u = "cm"
        if good_u is None and quantity == "time":
            good_u = "s"
        if return_quantity:
            unit_index = unit_list.index(good_u)
            # This avoids indexing errors
            if unit_index == 0:
                return self.quan(1, unit_list[0])
            # Number of orders of magnitude between unit and next one up
            OOMs = np.ceil(
                np.log10(
                    self.quan(1, unit_list[unit_index - 1])
                    / self.quan(1, unit_list[unit_index])
                )
            )
            # Backwards order of coefficients (e.g. [100, 10, 1])
            coeffs = 10 ** np.arange(OOMs)[::-1]
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
                self, dataset_type=self.dataset_type
            )
            # Now we do things that we need an instantiated index for
            # ...first off, we create our field_info now.
            oldsettings = np.geterr()
            np.seterr(all="ignore")
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
        for a in [
            "current_time",
            "domain_dimensions",
            "domain_left_edge",
            "domain_right_edge",
            "cosmological_simulation",
        ]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)
        if hasattr(self, "cosmological_simulation") and self.cosmological_simulation:
            for a in [
                "current_redshift",
                "omega_lambda",
                "omega_matter",
                "omega_radiation",
                "hubble_constant",
            ]:
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
        if "nbody" not in self.particle_types:
            mylog.debug("Creating Particle Union 'nbody'")
            ptypes = list(self.particle_types_raw)
            if hasattr(self, "_sph_ptypes"):
                for sph_ptype in self._sph_ptypes:
                    if sph_ptype in ptypes:
                        ptypes.remove(sph_ptype)
            if ptypes:
                nbody_ptypes = []
                for ptype in ptypes:
                    if (ptype, "particle_mass") in self.field_info:
                        nbody_ptypes.append(ptype)
                pu = ParticleUnion("nbody", nbody_ptypes)
                nfields = self.add_particle_union(pu)
                if nfields == 0:
                    mylog.debug("zero common fields, skipping particle union 'nbody'")
        self.field_info.setup_extra_union_fields()
        mylog.debug("Loading field plugins.")
        self.field_info.load_all_plugins(self.default_fluid_type)
        deps, unloaded = self.field_info.check_derived_fields()
        self.field_dependencies.update(deps)
        self.fields = FieldTypeContainer(self)
        self.index.field_list = sorted(self.field_list)
        self._last_freq = (None, None)

    def set_field_label_format(self, format_property, value):
        """
        Set format properties for how fields will be written
        out. Accepts

        format_property : string indicating what property to set
        value: the value to set for that format_property
        """
        available_formats = {"ionization_label": ("plus_minus", "roman_numeral")}
        if format_property in available_formats:
            if value in available_formats[format_property]:
                setattr(self, f"_{format_property}_format", value)
            else:
                raise ValueError(
                    "{0} not an acceptable value for format_property "
                    "{1}. Choices are {2}.".format(
                        value, format_property, available_formats[format_property]
                    )
                )
        else:
            raise ValueError(
                "{0} not a recognized format_property. Available"
                "properties are: {1}".format(
                    format_property, list(available_formats.keys())
                )
            )

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
            kwargs["ordering"] = ordering
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
            self.no_cgs_equiv_length = True
        elif self.geometry == "geographic":
            cls = GeographicCoordinateHandler
            self.no_cgs_equiv_length = True
        elif self.geometry == "internal_geographic":
            cls = InternalGeographicCoordinateHandler
            self.no_cgs_equiv_length = True
        elif self.geometry == "spectral_cube":
            cls = SpectralCubeCoordinateHandler
        else:
            raise YTGeometryNotSupported(self.geometry)
        self.coordinates = cls(self, **kwargs)

    def add_particle_union(self, union):
        # No string lookups here, we need an actual union.
        f = self.particle_fields_by_type

        # find fields common to all particle types in the union
        fields = set_intersection([f[s] for s in union if s in self.particle_types_raw])

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
        fields = [(union.name, field) for field in fields]
        new_fields = [_ for _ in fields if _ not in self.field_list]
        self.field_list.extend(new_fields)
        new_field_info_fields = [
            _ for _ in fields if _ not in self.field_info.field_list
        ]
        self.field_info.field_list.extend(new_field_info_fields)
        self.index.field_list = sorted(self.field_list)
        # Give ourselves a chance to add them here, first, then...
        # ...if we can't find them, we set them up as defaults.
        new_fields = self._setup_particle_types([union.name])
        self.field_info.find_dependencies(new_fields)
        return len(new_fields)

    def add_particle_filter(self, filter):
        """Add particle filter to the dataset.

        Add ``filter`` to the dataset and set up relavent derived_field.
        It will also add any ``filtered_type`` that the ``filter`` depends on.

        """
        # This requires an index
        self.index
        # This is a dummy, which we set up to enable passthrough of "all"
        # concatenation fields.
        n = getattr(filter, "name", filter)
        self.known_filters[n] = None
        if isinstance(filter, str):
            used = False
            f = filter_registry.get(filter, None)
            if f is None:
                return False
            used = self._setup_filtered_type(f)
            if used:
                filter = f
        else:
            used = self._setup_filtered_type(filter)
        if not used:
            self.known_filters.pop(n, None)
            return False
        self.known_filters[filter.name] = filter
        return True

    def _setup_filtered_type(self, filter):
        # Check if the filtered_type of this filter is known,
        # otherwise add it first if it is in the filter_registry
        if filter.filtered_type not in self.known_filters.keys():
            if filter.filtered_type in filter_registry:
                add_success = self.add_particle_filter(filter.filtered_type)
                if add_success:
                    mylog.info(
                        "Added filter dependency '%s' for '%s'",
                        filter.filtered_type,
                        filter.name,
                    )

        if not filter.available(self.derived_field_list):
            raise YTIllDefinedParticleFilter(
                filter, filter.missing(self.derived_field_list)
            )
        fi = self.field_info
        fd = self.field_dependencies
        available = False
        for fn in self.derived_field_list:
            if fn[0] == filter.filtered_type:
                # Now we can add this
                available = True
                self.derived_field_list.append((filter.name, fn[1]))
                fi[filter.name, fn[1]] = filter.wrap_func(fn, fi[fn])
                # Now we append the dependencies
                fd[filter.name, fn[1]] = fd[fn]
        if available:
            if filter.name not in self.particle_types:
                self.particle_types += (filter.name,)
            if filter.name not in self.filtered_particle_types:
                self.filtered_particle_types.append(filter.name)
            if hasattr(self, "_sph_ptypes"):
                if filter.filtered_type == self._sph_ptypes[0]:
                    mylog.warning(
                        "It appears that you are filtering on an SPH field "
                        "type. It is recommended to use 'gas' as the "
                        "filtered particle type in this case instead."
                    )
                if filter.filtered_type in (self._sph_ptypes + ("gas",)):
                    self._sph_ptypes = self._sph_ptypes + (filter.name,)
            new_fields = self._setup_particle_types([filter.name])
            deps, _ = self.field_info.check_derived_fields(new_fields)
            self.field_dependencies.update(deps)
        return available

    def _setup_particle_types(self, ptypes=None):
        df = []
        if ptypes is None:
            ptypes = self.ds.particle_types_raw
        for ptype in set(ptypes):
            df += self._setup_particle_type(ptype)
        return df

    _last_freq = (None, None)
    _last_finfo = None

    def _get_field_info(self, ftype, fname=None):
        self.index

        # store the original inputs in case we need to raise an error
        INPUT = ftype, fname
        if fname is None:
            try:
                ftype, fname = ftype.name
            except AttributeError:
                ftype, fname = "unknown", ftype

        # storing this condition before altering it
        guessing_type = ftype == "unknown"
        if guessing_type:
            ftype = self._last_freq[0] or ftype
        field = (ftype, fname)

        if (
            field == self._last_freq
            and field not in self.field_info.field_aliases.values()
        ):
            return self._last_finfo
        if field in self.field_info:
            self._last_freq = field
            self._last_finfo = self.field_info[(ftype, fname)]
            return self._last_finfo

        try:
            # Sometimes, if guessing_type == True, this will be switched for
            # the type of field it is.  So we look at the field type and
            # determine if we need to change the type.
            fi = self._last_finfo = self.field_info[fname]
            if (
                fi.sampling_type == "particle"
                and self._last_freq[0] not in self.particle_types
            ):
                field = "all", field[1]
            elif (
                not fi.sampling_type == "particle"
                and self._last_freq[0] not in self.fluid_types
            ):
                field = self.default_fluid_type, field[1]
            self._last_freq = field
            return self._last_finfo
        except KeyError:
            pass

        # We also should check "all" for particles, which can show up if you're
        # mixing deposition/gas fields with particle fields.
        if guessing_type:
            if hasattr(self, "_sph_ptype"):
                to_guess = [self.default_fluid_type, "all"]
            else:
                to_guess = ["all", self.default_fluid_type]
            to_guess += list(self.fluid_types) + list(self.particle_types)
            for ftype in to_guess:
                if (ftype, fname) in self.field_info:
                    self._last_freq = (ftype, fname)
                    self._last_finfo = self.field_info[(ftype, fname)]
                    return self._last_finfo
        raise YTFieldNotFound(field=INPUT, ds=self)

    def _setup_classes(self):
        # Called by subclass
        self.object_types = []
        self.objects = []
        self.plots = []
        for name, cls in sorted(data_object_registry.items()):
            if name in self._index_class._unsupported_objects:
                setattr(self, name, _unsupported_object(self, name))
                continue
            self._add_object_class(name, cls)
        self.object_types.sort()

    def _add_object_class(self, name, base):
        # skip projection data objects that don't make sense
        # for this type of data
        if "proj" in name and name != self._proj_type:
            return
        elif "proj" in name:
            name = "proj"
        self.object_types.append(name)
        obj = functools.partial(base, ds=weakref.proxy(self))
        obj.__doc__ = base.__doc__
        setattr(self, name, obj)

    def _find_extremum(self, field, ext, source=None, to_array=True):
        """
        Find the extremum value of a field in a data object (source) and its position.

        Parameters
        ----------
        field: str or tuple(str, str)

        ext: str
            'min' or 'max', select an extremum

        source: a Yt data object

        to_array: bool
            select the return type.

        Returns
        -------
        val, coords

        val: unyt.unyt_quantity
            extremum value detected

        coords: unyt.unyt_array or list(unyt.unyt_quantity)
            Conversion to a single unyt_array object is only possible for coordinate
            systems with homogeneous dimensions across axes (i.e. cartesian).
        """
        ext = ext.lower()
        if source is None:
            source = self.all_data()
        method = {
            "min": source.quantities.min_location,
            "max": source.quantities.max_location,
        }[ext]
        val, x1, x2, x3 = method(field)
        coords = [x1, x2, x3]
        mylog.info("%s value is %0.5e at %0.16f %0.16f %0.16f", ext, val, *coords)
        if to_array:
            if any(x.units.is_dimensionless for x in coords):
                mylog.warning(
                    f"dataset {self} has angular coordinates. "
                    "Use 'to_array=False' to preserve "
                    "dimensionality in each coordinate."
                )

            # force conversion to length
            alt_coords = []
            for x in coords:
                alt_coords.append(
                    self.quan(x.v, "code_length")
                    if x.units.is_dimensionless
                    else x.to("code_length")
                )
            coords = self.arr(alt_coords, dtype="float64").to("code_length")
        return val, coords

    def find_max(self, field, source=None, to_array=True):
        """
        Returns (value, location) of the maximum of a given field.

        This is a wrapper around _find_extremum
        """
        mylog.debug("Searching for maximum value of %s", field)
        return self._find_extremum(field, "max", source=source, to_array=to_array)

    def find_min(self, field, source=None, to_array=True):
        """
        Returns (value, location) for the minimum of a given field.

        This is a wrapper around _find_extremum
        """
        mylog.debug("Searching for minimum value of %s", field)
        return self._find_extremum(field, "min", source=source, to_array=to_array)

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
            funit = self._get_field_info(field).units
            out.append(self.arr(np.empty((len(coords),)), funit))
            for coord_index, coord in enumerate(coords):
                out[field_index][coord_index] = self.point(coord)[field]
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
        self.index
        if find_max:
            c = self.find_max("density")[1]
        else:
            c = (self.domain_right_edge + self.domain_left_edge) / 2.0
        return self.region(c, self.domain_left_edge, self.domain_right_edge, **kwargs)

    def box(self, left_edge, right_edge, **kwargs):
        """
        box is a wrapper to the Region object for creating a region
        without having to specify a *center* value.  It assumes the center
        is the midpoint between the left_edge and right_edge.

        Keyword arguments are passed to the initializer of the YTRegion object
        (e.g. ds.region).
        """
        # we handle units in the region data object
        # but need to check if left_edge or right_edge is a
        # list or other non-array iterable before calculating
        # the center
        if isinstance(left_edge[0], YTQuantity):
            left_edge = YTArray(left_edge)
            right_edge = YTArray(right_edge)

        left_edge = np.asanyarray(left_edge, dtype="float64")
        right_edge = np.asanyarray(right_edge, dtype="float64")
        c = (left_edge + right_edge) / 2.0
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
        if not self.particles_exist:
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
        return self.refine_by ** (l1 - l0)

    def _assign_unit_system(self, unit_system):
        if unit_system == "cgs":
            current_mks_unit = None
        else:
            current_mks_unit = "A"
        magnetic_unit = getattr(self, "magnetic_unit", None)
        if magnetic_unit is not None:
            if unit_system == "mks":
                if current_mks not in self.magnetic_unit.units.dimensions.free_symbols:
                    self.magnetic_unit = self.magnetic_unit.to("gauss").to("T")
                self.unit_registry.modify("code_magnetic", self.magnetic_unit.value)
            else:
                # if the magnetic unit is in T, we need to create the code unit
                # system as an MKS-like system
                if current_mks in self.magnetic_unit.units.dimensions.free_symbols:
                    self.magnetic_unit = self.magnetic_unit.to("T").to("gauss")
                # The following modification ensures that we get the conversion to
                # cgs correct
                self.unit_registry.modify(
                    "code_magnetic", self.magnetic_unit.value * 0.1 ** 0.5
                )

        us = create_code_unit_system(
            self.unit_registry, current_mks_unit=current_mks_unit
        )
        if unit_system != "code":
            us = unit_system_registry[str(unit_system).lower()]
        self.unit_system = us
        self.unit_registry.unit_system = self.unit_system

    def _create_unit_registry(self, unit_system):
        from yt.units import dimensions

        # yt assumes a CGS unit system by default (for back compat reasons).
        # Since unyt is MKS by default we specify the MKS values of the base
        # units in the CGS system. So, for length, 1 cm = .01 m. And so on.
        self.unit_registry = UnitRegistry(unit_system=unit_system)
        self.unit_registry.add("code_length", 0.01, dimensions.length)
        self.unit_registry.add("code_mass", 0.001, dimensions.mass)
        self.unit_registry.add("code_density", 1000.0, dimensions.density)
        self.unit_registry.add(
            "code_specific_energy", 1.0, dimensions.energy / dimensions.mass
        )
        self.unit_registry.add("code_time", 1.0, dimensions.time)
        if unit_system == "mks":
            self.unit_registry.add("code_magnetic", 1.0, dimensions.magnetic_field)
        else:
            self.unit_registry.add(
                "code_magnetic", 0.1 ** 0.5, dimensions.magnetic_field_cgs
            )
        self.unit_registry.add("code_temperature", 1.0, dimensions.temperature)
        self.unit_registry.add("code_pressure", 0.1, dimensions.pressure)
        self.unit_registry.add("code_velocity", 0.01, dimensions.velocity)
        self.unit_registry.add("code_metallicity", 1.0, dimensions.dimensionless)
        self.unit_registry.add("h", 1.0, dimensions.dimensionless, r"h")
        self.unit_registry.add("a", 1.0, dimensions.dimensionless)

    def set_units(self):
        """
        Creates the unit registry for this dataset.

        """

        if getattr(self, "cosmological_simulation", False):
            # this dataset is cosmological, so add cosmological units.
            self.unit_registry.modify("h", self.hubble_constant)
            if getattr(self, "current_redshift", None) is not None:
                # Comoving lengths
                for my_unit in ["m", "pc", "AU", "au"]:
                    new_unit = f"{my_unit}cm"
                    my_u = Unit(my_unit, registry=self.unit_registry)
                    self.unit_registry.add(
                        new_unit,
                        my_u.base_value / (1 + self.current_redshift),
                        dimensions.length,
                        "\\rm{%s}/(1+z)" % my_unit,
                        prefixable=True,
                    )
                self.unit_registry.modify("a", 1 / (1 + self.current_redshift))

        self.set_code_units()

    def setup_cosmology(self):
        """
        If this dataset is cosmological, add a cosmology object.
        """
        if not getattr(self, "cosmological_simulation", False):
            return

        # Set dynamical dark energy parameters
        use_dark_factor = getattr(self, "use_dark_factor", False)
        w_0 = getattr(self, "w_0", -1.0)
        w_a = getattr(self, "w_a", 0.0)

        # many frontends do not set this
        setdefaultattr(self, "omega_radiation", 0.0)

        self.cosmology = Cosmology(
            hubble_constant=self.hubble_constant,
            omega_matter=self.omega_matter,
            omega_lambda=self.omega_lambda,
            omega_radiation=self.omega_radiation,
            use_dark_factor=use_dark_factor,
            w_0=w_0,
            w_a=w_a,
        )

        if getattr(self, "current_redshift", None) is not None:
            self.critical_density = self.cosmology.critical_density(
                self.current_redshift
            )
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
        vel_unit = getattr(self, "velocity_unit", self.length_unit / self.time_unit)
        pressure_unit = getattr(
            self,
            "pressure_unit",
            self.mass_unit / (self.length_unit * (self.time_unit) ** 2),
        )
        temperature_unit = getattr(self, "temperature_unit", 1.0)
        density_unit = getattr(
            self, "density_unit", self.mass_unit / self.length_unit ** 3
        )
        specific_energy_unit = getattr(self, "specific_energy_unit", vel_unit ** 2)
        self.unit_registry.modify("code_velocity", vel_unit)
        self.unit_registry.modify("code_temperature", temperature_unit)
        self.unit_registry.modify("code_pressure", pressure_unit)
        self.unit_registry.modify("code_density", density_unit)
        self.unit_registry.modify("code_specific_energy", specific_energy_unit)
        # domain_width does not yet exist
        if self.domain_left_edge is not None and self.domain_right_edge is not None:
            DW = self.arr(self.domain_right_edge - self.domain_left_edge, "code_length")
            self.unit_registry.add(
                "unitary", float(DW.max() * DW.units.base_value), DW.units.dimensions
            )

    @classmethod
    def _validate_units_override_keys(cls, units_override):
        valid_keys = set(cls.default_units.keys())
        invalid_keys_found = set(units_override.keys()) - valid_keys
        if invalid_keys_found:
            raise ValueError(
                f"units_override contains invalid keys: {invalid_keys_found}"
            )

    default_units = {
        "length_unit": "cm",
        "time_unit": "s",
        "mass_unit": "g",
        "velocity_unit": "cm/s",
        "magnetic_unit": "gauss",
        "temperature_unit": "K",
    }

    @classmethod
    def _sanitize_units_override(cls, units_override):
        """
        Convert units_override values to valid input types for unyt.
        Throw meaningful errors early if units_override is ill-formed.

        Parameters
        ----------
        units_override : dict

            keys should be strings with format  "<dim>_unit" (e.g. "mass_unit"), and
            need to match a key in cls.default_units

            values should be mappable to unyt.unyt_quantity objects, and can be any
            combinations of:
                - unyt.unyt_quantity
                - 2-long sequence (tuples, list, ...) with types (number, str)
                  e.g. (10, "km"), (0.1, "s")
                - number (in which case the associated is taken from cls.default_unit)


        Raises
        ------
        TypeError
            If unit_override has invalid types

        ValueError
            If provided units do not match the intended dimensionality,
            or in case of a zero scaling factor.

        """
        uo = {}
        if units_override is None:
            return uo

        cls._validate_units_override_keys(units_override)

        for key in cls.default_units:
            try:
                val = units_override[key]
            except KeyError:
                continue

            # Now attempt to instanciate a unyt.unyt_quantity from val ...
            try:
                # ... directly (valid if val is a number, or a unyt_quantity)
                uo[key] = YTQuantity(val)
                continue
            except RuntimeError:
                # note that unyt.unyt_quantity throws RuntimeError in lieu of TypeError
                pass
            try:
                # ... with tuple unpacking (valid if val is a sequence)
                uo[key] = YTQuantity(*val)
                continue
            except (RuntimeError, TypeError, UnitParseError):
                pass
            raise TypeError(
                "units_override values should be 2-sequence (float, str), "
                "YTQuantity objects or real numbers; "
                f"received {val} with type {type(val)}."
            )
        for key, q in uo.items():
            if q.units.is_dimensionless:
                uo[key] = YTQuantity(q, cls.default_units[key])
            try:
                uo[key].to(cls.default_units[key])
            except UnitConversionError as err:
                raise ValueError(
                    "Inconsistent dimensionality in units_override. "
                    f"Received {key} = {uo[key]}"
                ) from err
            if 1 / uo[key].value == np.inf:
                raise ValueError(
                    f"Invalid 0 normalisation factor in units_override for {key}."
                )
        return uo

    def _override_code_units(self):
        if not self.units_override:
            return

        mylog.warning(
            "Overriding code units: Use this option only if you know that the "
            "dataset doesn't define the units correctly or at all."
        )
        for ukey, val in self.units_override.items():
            mylog.info("Overriding %s: %s.", ukey, val)
            setattr(self, ukey, self.quan(val))

    _units = None
    _unit_system_id = None

    @property
    def units(self):
        current_uid = self.unit_registry.unit_system_id
        if self._units is not None and self._unit_system_id == current_uid:
            return self._units
        self._unit_system_id = current_uid
        self._units = UnitContainer(self.unit_registry)
        return self._units

    _arr = None

    @property
    def arr(self):
        """Converts an array into a :class:`yt.units.yt_array.YTArray`

        The returned YTArray will be dimensionless by default, but can be
        cast to arbitrary units using the ``units`` keyword argument.

        Parameters
        ----------

        input_array : Iterable
            A tuple, list, or array to attach units to
        units: String unit specification, unit symbol or astropy object
            The units of the array. Powers must be specified using python syntax
            (cm**3, not cm^3).
        input_units : Deprecated in favor of 'units'
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
        self._arr = functools.partial(YTArray, registry=self.unit_registry)
        return self._arr

    _quan = None

    @property
    def quan(self):
        """Converts an scalar into a :class:`yt.units.yt_array.YTQuantity`

        The returned YTQuantity will be dimensionless by default, but can be
        cast to arbitrary units using the ``units`` keyword argument.

        Parameters
        ----------

        input_scalar : an integer or floating point scalar
            The scalar to attach units to
        units: String unit specification, unit symbol or astropy object
            The units of the quantity. Powers must be specified using python
            syntax (cm**3, not cm^3).
        input_units : Deprecated in favor of 'units'
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

    def add_field(self, name, function, sampling_type, **kwargs):
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
        sampling_type: str
           "cell" or "particle" or "local"
        units : str
           A plain text string encoding the unit.  Powers must be in
           python syntax (** instead of ^).
        take_log : bool
           Describes whether the field should be logged
        validators : list
           A list of :class:`FieldValidator` objects
        vector_field : bool
           Describes the dimensionality of the field.  Currently unused.
        display_name : str
           A name used in the plots
        force_override : bool
           Whether to override an existing derived field. Does not work with
           on-disk fields.

        """
        self.index
        override = kwargs.get("force_override", False)
        if override and name in self.index.field_list:
            raise RuntimeError(
                "force_override is only meant to be used with "
                "derived fields, not on-disk fields."
            )
        # Handle the case where the field has already been added.
        if not override and name in self.field_info:
            mylog.warning(
                "Field %s already exists. To override use `force_override=True`.", name,
            )

        self.field_info.add_field(name, function, sampling_type, **kwargs)
        self.field_info._show_field_errors.append(name)
        deps, _ = self.field_info.check_derived_fields([name])
        self.field_dependencies.update(deps)

    def add_mesh_sampling_particle_field(self, sample_field, ptype="all"):
        """Add a new mesh sampling particle field

        Creates a new particle field which has the value of the
        *deposit_field* at the location of each particle of type
        *ptype*.

        Parameters
        ----------

        sample_field : tuple
           The field name tuple of the mesh field to be deposited onto
           the particles. This must be a field name tuple so yt can
           appropriately infer the correct particle type.
        ptype : string, default 'all'
           The particle type onto which the deposition will occur.

        Returns
        -------

        The field name tuple for the newly created field.

        Examples
        --------
        >>> ds = yt.load('output_00080/info_00080.txt')
        ... ds.add_mesh_sampling_particle_field(('gas', 'density'), ptype='all')

        >>> print('The density at the location of the particle is:')
        ... print(ds.r['all', 'cell_gas_density'])
        The density at the location of the particle is:
        [9.33886124e-30 1.22174333e-28 1.20402333e-28 ... 2.77410331e-30
         8.79467609e-31 3.50665136e-30] g/cm**3

        >>> len(ds.r['all', 'cell_gas_density']) == len(ds.r['all', 'particle_ones'])
        True

        """
        if isinstance(sample_field, tuple):
            ftype, sample_field = sample_field[0], sample_field[1]
        else:
            raise RuntimeError

        return self.index._add_mesh_sampling_particle_field(sample_field, ftype, ptype)

    def add_deposited_particle_field(
        self, deposit_field, method, kernel_name="cubic", weight_field="particle_mass"
    ):
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
           `nearest` and `count`.
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

        units = self.field_info[ptype, deposit_field].output_units
        take_log = self.field_info[ptype, deposit_field].take_log
        name_map = {
            "sum": "sum",
            "std": "std",
            "cic": "cic",
            "weighted_mean": "avg",
            "nearest": "nn",
            "simple_smooth": "ss",
            "count": "count",
        }
        field_name = "%s_" + name_map[method] + "_%s"
        field_name = field_name % (ptype, deposit_field.replace("particle_", ""))

        if method == "count":
            field_name = f"{ptype}_count"
            if ("deposit", field_name) in self.field_info:
                mylog.warning("The deposited field %s already exists", field_name)
                return ("deposit", field_name)
            else:
                units = "dimensionless"
                take_log = False

        def _deposit_field(field, data):
            """
            Create a grid field for particle quantities using given method.
            """
            pos = data[ptype, "particle_position"]
            fields = [data[ptype, deposit_field]]
            if method == "weighted_mean":
                fields.append(data[ptype, weight_field])
            fields = [np.ascontiguousarray(f) for f in fields]
            d = data.deposit(pos, fields, method=method, kernel_name=kernel_name)
            d = data.ds.arr(d, units=units)
            if method == "weighted_mean":
                d[np.isnan(d)] = 0.0
            return d

        self.add_field(
            ("deposit", field_name),
            function=_deposit_field,
            sampling_type="cell",
            units=units,
            take_log=take_log,
            validators=[ValidateSpatial()],
        )
        return ("deposit", field_name)

    def add_smoothed_particle_field(
        self, smooth_field, method="volume_weighted", nneighbors=64, kernel_name="cubic"
    ):
        """Add a new smoothed particle field

        WARNING: This method is deprecated since yt-4.0.

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
        kernel_name : string, default `cubic`
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        issue_deprecation_warning("This method is deprecated. " + DEP_MSG_SMOOTH_FIELD)

    def add_gradient_fields(self, fields=None, input_field=None):
        """Add gradient fields.

        Creates four new grid-based fields that represent the components of the gradient
        of an existing field, plus an extra field for the magnitude of the gradient. The
        gradient is computed using second-order centered differences.

        Parameters
        ----------
        fields : str or tuple(str, str), or a list of the previous
            Label(s) for at least one field. Can either represent a tuple
            (<field type>, <field fname>) or simply the field name.
            Warning: several field types may match the provided field name,
            in which case the first one discovered internally is used.

        Returns
        -------
        A list of field name tuples for the newly created fields.

        Raises
        ------
        YTFieldNotParsable
            If fields are not parsable to yt field keys.

        YTFieldNotFound :
            If at least one field can not be identified.

        Examples
        --------

        >>> grad_fields = ds.add_gradient_fields(("gas","density"))
        >>> print(grad_fields)
        ... [('gas', 'density_gradient_x'), ('gas', 'density_gradient_y'),
        ...  ('gas', 'density_gradient_z'), ('gas', 'density_gradient_magnitude')]

        Note that the above example assumes ds.geometry == 'cartesian'. In general,
        the function will create gradient components along the axes of the dataset
        coordinate system.
        For instance, with cylindrical data, one gets 'density_gradient_<r,theta,z>'

        """
        if input_field is not None:
            issue_deprecation_warning(
                "keyword argument 'input_field' is deprecated in favor of 'fields' "
                "and will be removed in a future version of yt."
            )
            if fields is not None:
                raise TypeError(
                    "Can not use both 'fields' and 'input_field' keyword arguments"
                )
            fields = input_field
        if fields is None:
            raise TypeError("Missing required positional argument: fields")

        self.index
        data_obj = self.all_data()
        explicit_fields = data_obj._determine_fields(fields)
        grad_fields = []
        for ftype, fname in explicit_fields:
            units = self.field_info[ftype, fname].units
            setup_gradient_fields(self.field_info, (ftype, fname), units)
            # Now we make a list of the fields that were just made, to check them
            # and to return them
            grad_fields += [
                (ftype, f"{fname}_gradient_{suffix}")
                for suffix in self.coordinates.axis_order
            ]
            grad_fields.append((ftype, f"{fname}_gradient_magnitude"))
            deps, _ = self.field_info.check_derived_fields(grad_fields)
            self.field_dependencies.update(deps)
        return grad_fields

    _max_level = None

    @property
    def max_level(self):
        if self._max_level is None:
            self._max_level = self.index.max_level
        return self._max_level

    @max_level.setter
    def max_level(self, value):
        self._max_level = value

    _min_level = None

    @property
    def min_level(self):
        if self._min_level is None:
            self._min_level = self.index.min_level
        return self._min_level

    @min_level.setter
    def min_level(self, value):
        self._min_level = value

    def define_unit(self, symbol, value, tex_repr=None, offset=None, prefixable=False):
        """
        Define a new unit and add it to the dataset's unit registry.

        Parameters
        ----------
        symbol : string
            The symbol for the new unit.
        value : tuple or ~yt.units.yt_array.YTQuantity
            The definition of the new unit in terms of some other units. For example,
            one would define a new "mph" unit with (1.0, "mile/hr")
        tex_repr : string, optional
            The LaTeX representation of the new unit. If one is not supplied, it will
            be generated automatically based on the symbol string.
        offset : float, optional
            The default offset for the unit. If not set, an offset of 0 is assumed.
        prefixable : bool, optional
            Whether or not the new unit can use SI prefixes. Default: False

        Examples
        --------
        >>> ds.define_unit("mph", (1.0, "mile/hr"))
        >>> two_weeks = YTQuantity(14.0, "days")
        >>> ds.define_unit("fortnight", two_weeks)
        """
        define_unit(
            symbol,
            value,
            tex_repr=tex_repr,
            offset=offset,
            prefixable=prefixable,
            registry=self.unit_registry,
        )


def _reconstruct_ds(*args, **kwargs):
    datasets = ParameterFileStore()
    ds = datasets.get_ds_hash(*args)
    return ds


@functools.total_ordering
class ParticleFile:
    def __init__(self, ds, io, filename, file_id, range=None):
        self.ds = ds
        self.io = weakref.proxy(io)
        self.filename = filename
        self.file_id = file_id
        if range is None:
            range = (None, None)
        self.start, self.end = range
        self.total_particles = self.io._count_particles(self)
        # Now we adjust our start/end, in case there are fewer particles than
        # we realized
        if self.start is None:
            self.start = 0
        self.end = max(self.total_particles.values()) + self.start

    def select(self, selector):
        pass

    def count(self, selector):
        pass

    def _calculate_offsets(self, fields, pcounts):
        pass

    def __lt__(self, other):
        if self.filename != other.filename:
            return self.filename < other.filename
        return self.start < other.start

    def __eq__(self, other):
        if self.filename != other.filename:
            return False
        return self.start == other.start

    def __hash__(self):
        return hash((self.filename, self.file_id, self.start, self.end))


class ParticleDataset(Dataset):
    _unit_base = None
    filter_bbox = False
    _proj_type = "particle_proj"

    def __init__(
        self,
        filename,
        dataset_type=None,
        file_style=None,
        units_override=None,
        unit_system="cgs",
        index_order=None,
        index_filename=None,
    ):
        self.index_order = validate_index_order(index_order)
        self.index_filename = index_filename
        super(ParticleDataset, self).__init__(
            filename,
            dataset_type=dataset_type,
            file_style=file_style,
            units_override=units_override,
            unit_system=unit_system,
        )


def validate_index_order(index_order):
    if index_order is None:
        index_order = (7, 5)
    elif not iterable(index_order):
        index_order = (int(index_order), 1)
    else:
        if len(index_order) != 2:
            raise RuntimeError(
                "Tried to load a dataset with index_order={}, but "
                "index_order\nmust be an integer or a two-element tuple of "
                "integers.".format(index_order)
            )
        index_order = tuple([int(o) for o in index_order])
    return index_order

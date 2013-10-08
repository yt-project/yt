"""
Generalized Enzo output objects, both static and time-series.

Presumably at some point EnzoRun will be absorbed into here.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import string, re, gc, time, os, os.path, weakref

from yt.funcs import *

from yt.config import ytcfg
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.utilities.parameter_file_storage import \
    ParameterFileStore, \
    NoParameterShelf, \
    output_type_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.data_objects.particle_filters import \
    filter_registry
from yt.data_objects.particle_unions import \
    ParticleUnion
from yt.utilities.minimal_representation import \
    MinimalStaticOutput

from yt.geometry.coordinate_handler import \
    CartesianCoordinateHandler, \
    PolarCoordinateHandler, \
    CylindricalCoordinateHandler

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is contant up
# to here, and then have it instantiate EnzoStaticOutputs as appropriate.

_cached_pfs = weakref.WeakValueDictionary()
_pf_store = ParameterFileStore()

class StaticOutput(object):

    default_fluid_type = "gas"
    fluid_types = ("gas","deposit")
    particle_types = ("io",) # By default we have an 'all'
    particle_types_raw = ("io",)
    geometry = "cartesian"
    coordinates = None
    max_level = 99
    storage_filename = None
    _particle_mass_name = None
    _particle_coordinates_name = None
    _particle_velocity_name = None
    particle_unions = None
    known_filters = None

    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            output_type_registry[name] = cls
            mylog.debug("Registering: %s as %s", name, cls)

    def __new__(cls, filename=None, *args, **kwargs):
        if not isinstance(filename, types.StringTypes):
            obj = object.__new__(cls)
            # The Stream frontend uses a StreamHandler object to pass metadata
            # to __init__.
            is_stream = (hasattr(filename, 'get_fields') and
                         hasattr(filename, 'get_particle_type'))
            if not is_stream:
                obj.__init__(filename, *args, **kwargs)
            return obj
        apath = os.path.abspath(filename)
        if not os.path.exists(apath): raise IOError(filename)
        if apath not in _cached_pfs:
            obj = object.__new__(cls)
            if obj._skip_cache is False:
                _cached_pfs[apath] = obj
        else:
            obj = _cached_pfs[apath]
        return obj

    def __init__(self, filename, data_style=None, file_style=None):
        """
        Base class for generating new output types.  Principally consists of
        a *filename* and a *data_style* which will be passed on to children.
        """
        self.data_style = data_style
        self.file_style = file_style
        self.conversion_factors = {}
        self.parameters = {}
        self.known_filters = self.known_filters or {}
        self.particle_unions = self.particle_unions or {}

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

        self._parse_parameter_file()
        self._setup_coordinate_handler()
        self._set_units()
        self._set_derived_attrs()

        # Because we need an instantiated class to check the pf's existence in
        # the cache, we move that check to here from __new__.  This avoids
        # double-instantiation.
        try:
            _pf_store.check_pf(self)
        except NoParameterShelf:
            pass
        self.print_key_parameters()

        self.create_field_info()

    def _set_derived_attrs(self):
        self.domain_center = 0.5 * (self.domain_right_edge + self.domain_left_edge)
        self.domain_width = self.domain_right_edge - self.domain_left_edge

    def __reduce__(self):
        args = (self._hash(),)
        return (_reconstruct_pf, args)

    def __repr__(self):
        return self.basename

    def _hash(self):
        s = "%s;%s;%s" % (self.basename,
            self.current_time, self.unique_identifier)
        try:
            import hashlib
            return hashlib.md5(s).hexdigest()
        except ImportError:
            return s.replace(";", "*")

    @property
    def _mrep(self):
        return MinimalStaticOutput(self)

    @property
    def _skip_cache(self):
        return False

    def hub_upload(self):
        self._mrep.upload()

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return False

    def __getitem__(self, key):
        """ Returns units, parameters, or conversion_factors in that order. """
        for d in [self.units, self.time_units, self.parameters, \
                  self.conversion_factors]:
            if key in d: return d[key]
        raise KeyError(key)

    def keys(self):
        """
        Returns a list of possible keys, from units, parameters and
        conversion_factors.

        """
        return self.units.keys() \
             + self.time_units.keys() \
             + self.parameters.keys() \
             + self.conversion_factors.keys()

    def __iter__(self):
        for ll in [self.units, self.time_units,
                   self.parameters, self.conversion_factors]:
            for i in ll.keys(): yield i

    def get_smallest_appropriate_unit(self, v):
        max_nu = 1e30
        good_u = None
        for unit in ['mpc', 'kpc', 'pc', 'au', 'rsun', 'km', 'cm']:
            vv = v*self[unit]
            if vv < max_nu and vv > 1.0:
                good_u = unit
                max_nu = v*self[unit]
        if good_u is None : good_u = 'cm'
        return good_u

    def has_key(self, key):
        """
        Checks units, parameters, and conversion factors. Returns a boolean.

        """
        return key in self.units or \
               key in self.time_units or \
               key in self.parameters or \
               key in self.conversion_factors

    _instantiated_hierarchy = None
    @property
    def hierarchy(self):
        if self._instantiated_hierarchy is None:
            if self._hierarchy_class == None:
                raise RuntimeError("You should not instantiate StaticOutput.")
            self._instantiated_hierarchy = self._hierarchy_class(
                self, data_style=self.data_style)
            # Now we do things that we need an instantiated hierarchy for
            if "all" not in self.particle_types:
                pu = ParticleUnion("all", list(self.particle_types_raw))
                self.add_particle_union(pu)
        return self._instantiated_hierarchy
    h = hierarchy  # alias

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

    def create_field_info(self):
        if getattr(self, "field_info", None) is None:
            # The setting up of fields occurs in the hierarchy, which is only
            # instantiated once.  So we have to double check to make sure that,
            # in the event of double-loads of a parameter file, we do not blow
            # away the exising field_info.
            self.field_info = FieldInfoContainer.create_with_fallback(
                                self._fieldinfo_fallback)
            self.field_info.update(self.coordinates.coordinate_fields())
        if getattr(self, "field_dependencies", None) is None:
            self.field_dependencies = {}

    def _setup_coordinate_handler(self):
        if self.geometry == "cartesian":
            self.coordinates = CartesianCoordinateHandler(self)
        elif self.geometry == "cylindrical":
            self.coordinates = CylindricalCoordinateHandler(self)
        elif self.geometry == "polar":
            self.coordinates = PolarCoordinateHandler(self)
        else:
            raise YTGeometryNotSupported(self.geometry)

    def add_particle_filter(self, filter):
        # This is a dummy, which we set up to enable passthrough of "all"
        # concatenation fields.
        n = getattr(filter, "name", filter)
        self.known_filters[n] = None
        if isinstance(filter, types.StringTypes):
            used = False
            for f in filter_registry[filter]:
                used = self.h._setup_filtered_type(f)
                if used:
                    filter = f
                    break
        else:
            used = self.h._setup_filtered_type(filter)
        if not used:
            self.known_filters.pop(n, None)
            return False
        self.known_filters[filter.name] = filter
        return True

    _last_freq = (None, None)
    _last_finfo = None
    def _get_field_info(self, ftype, fname):
        guessing_type = False
        if ftype == "unknown" and self._last_freq[0] != None:
            ftype = self._last_freq[0]
            guessing_type = True
        field = (ftype, fname)
        if field == self._last_freq:
            return self._last_finfo
        if field in self.field_info:
            self._last_freq = field
            self._last_finfo = self.field_info[(ftype, fname)]
            return self._last_finfo
        if fname == self._last_freq[1]:
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
        if guessing_type and ("all", fname) in self.field_info:
            self._last_freq = ("all", fname)
            self._last_finfo = self.field_info["all", fname]
            return self._last_finfo
        raise YTFieldNotFound((ftype, fname), self)

    def add_particle_union(self, union):
        # No string lookups here, we need an actual union.
        f = self.particle_fields_by_type
        fields = set_intersection([f[s] for s in union
                                   if s in self.particle_types_raw])
        self.particle_types += (union.name,)
        self.particle_unions[union.name] = union
        fields = [ (union.name, field) for field in fields]
        self.h.field_list.extend(fields)
        self.h._setup_unknown_fields(fields)
        self.h._setup_particle_types([union.name])

    def _setup_particle_type(self, ptype):
        mylog.debug("Don't know what to do with %s", ptype)
        return []

    @property
    def particle_fields_by_type(self):
        fields = defaultdict(list)
        for field in self.h.field_list:
            if field[0] in self.particle_types_raw:
                fields[field[0]].append(field[1])
        return fields

    @property
    def ires_factor(self):
        o2 = np.log2(self.refine_by)
        if o2 != int(o2):
            raise RuntimeError
        return int(o2)

def _reconstruct_pf(*args, **kwargs):
    pfs = ParameterFileStore()
    pf = pfs.get_pf_hash(*args)
    return pf


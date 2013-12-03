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

import string, re, gc, time, os, weakref, copy

from yt.funcs import *

from yt.config import ytcfg
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.utilities.parameter_file_storage import \
    ParameterFileStore, \
    NoParameterShelf, \
    output_type_registry
from yt.data_objects.data_containers import \
    data_object_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.data_objects.particle_filters import \
    filter_registry
from yt.data_objects.particle_unions import \
    ParticleUnion
from yt.utilities.minimal_representation import \
    MinimalDataset
from yt.utilities.io_handler import io_registry

from yt.geometry.coordinate_handler import \
    CartesianCoordinateHandler, \
    PolarCoordinateHandler, \
    CylindricalCoordinateHandler

# We want to support the movie format in the future.
# When such a thing comes to pass, I'll move all the stuff that is contant up
# to here, and then have it instantiate EnzoDatasets as appropriate.

_cached_pfs = weakref.WeakValueDictionary()
_pf_store = ParameterFileStore()

class Dataset(object):

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
    _unsupported_objects = ()
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

    def __init__(self, filename, dataset_type=None, file_style=None):
        """
        Base class for generating new output types.  Principally consists of
        a *filename* and a *dataset_type* which will be passed on to children.
        """
        self.dataset_type = dataset_type
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

        # Must be defined in subclass
        mylog.debug("Setting up classes.")
        self._setup_classes()

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
        return MinimalDataset(self)

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

    _instantiated_index = None
    @property
    def index(self):
        if self._instantiated_index is None:
            self._initialize_index()
        return self._instantiated_index

    def _initialize_index(self):
        index = self._index_class(
                self, dataset_type=self.dataset_type)
        self._instantiated_index = index

        # Here is where we initialize the index.
        mylog.debug("Initializing data IO")
        self._setup_data_io()

        mylog.debug("Detecting fields.")
        self._detect_fields()

        mylog.debug("Detecting fields in backup.")
        self._detect_fields_backup()

        mylog.debug("Adding unknown detected fields")
        self._setup_unknown_fields()

        mylog.debug("Setting up derived fields")
        self._setup_derived_fields()

        mylog.debug("Setting up particle fields")
        self._setup_particle_types()

        mylog.debug("Checking derived fields again")
        self._derived_fields_add(self._check_later)

        if "all" not in self.particle_types:
            mylog.debug("Creating Particle Union 'all'")
            pu = ParticleUnion("all", list(self.particle_types_raw))
            self.add_particle_union(pu)

    @property
    def hierarchy(self):
        return self.index

    @property
    def h(self):
        return self

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
    def _get_field_info(self, ftype, fname = None):
        self.index # This requires an index creation.
        if fname is None:
            ftype, fname = "unknown", ftype
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
        if guessing_type:
            for ftype in ("all", self.default_fluid_type):
                if (ftype, fname) in self.field_info:
                    self._last_freq = (ftype, fname)
                    self._last_finfo = self.field_info[(ftype, fname)]
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
        self.field_list.extend(fields)
        # Give ourselves a chance to add them here, first, then...
        # ...if we can't find them, we set them up as defaults.
        self._setup_particle_types([union.name])
        self._setup_unknown_fields(fields, self.field_info,
                                     skip_removal = True)

    def _setup_particle_type(self, ptype):
        mylog.debug("Don't know what to do with %s", ptype)
        return []

    @property
    def particle_fields_by_type(self):
        fields = defaultdict(list)
        for field in self.field_list:
            if field[0] in self.particle_types_raw:
                fields[field[0]].append(field[1])
        return fields

    @property
    def ires_factor(self):
        o2 = np.log2(self.refine_by)
        if o2 != int(o2):
            raise RuntimeError
        return int(o2)

    def relative_refinement(self, l0, l1):
        return self.refine_by**(l1-l0)

    def find_max(self, field):
        """
        Returns (value, center) of location of maximum for a given field.
        """
        mylog.debug("Searching for maximum value of %s", field)
        source = self.all_data()
        max_val, maxi, mx, my, mz = \
            source.quantities["MaxLocation"](field)
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f", 
              max_val, mx, my, mz)
        return max_val, np.array([mx, my, mz], dtype="float64")

    def _add_object_class(self, name, class_name, base, dd):
        self.object_types.append(name)
        obj = type(class_name, (base,), dd)
        setattr(self, name, obj)

    def _get_data_reader_dict(self):
        dd = { 'pf' : self }
        return dd

    def _setup_classes(self):
        # Called by subclass
        dd = self._get_data_reader_dict()
        self.object_types = []
        self.objects = []
        self.plots = []
        for name, cls in sorted(data_object_registry.items()):
            if name in self._unsupported_objects:
                setattr(self, name,
                    _unsupported_object(self.parameter_file, name))
                continue
            cname = cls.__name__
            if cname.endswith("Base"): cname = cname[:-4]
            self._add_object_class(name, cname, cls, dd)
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

    def _detect_fields_backup(self):
        # grab fields from backup file as well, if present
        return
        try:
            backup_filename = self.parameter_file.backup_filename
            f = h5py.File(backup_filename, 'r')
            g = f["data"]
            grid = self.grids[0] # simply check one of the grids
            grid_group = g["grid_%010i" % (grid.id - grid._id_offset)]
            for field_name in grid_group:
                if field_name != 'particles':
                    self.field_list.append(field_name)
        except KeyError:
            return
        except IOError:
            return

    def _setup_particle_types(self, ptypes = None):
        mname = self._particle_mass_name
        cname = self._particle_coordinates_name
        vname = self._particle_velocity_name
        # We require overriding if any of this is true
        df = []
        if ptypes is None: ptypes = self.particle_types_raw
        if None in (mname, cname, vname): 
            # If we don't know what to do, then let's not.
            for ptype in set(ptypes):
                df += self._setup_particle_type(ptype)
            # Now we have a bunch of new fields to add!
            # This is where the dependencies get calculated.
            self._derived_fields_add(df)
            return
        fi = self.field_info
        def _get_conv(cf):
            def _convert(data):
                return data.convert(cf)
            return _convert
        for ptype in ptypes:
            fi.add_field((ptype, vname), function=NullFunc,
                particle_type = True,
                convert_function=_get_conv("velocity"),
                units = r"\mathrm{cm}/\mathrm{s}")
            df.append((ptype, vname))
            fi.add_field((ptype, mname), function=NullFunc,
                particle_type = True,
                convert_function=_get_conv("mass"),
                units = r"\mathrm{g}")
            df.append((ptype, mname))
            df += particle_deposition_functions(ptype, cname, mname, fi)
            df += particle_scalar_functions(ptype, cname, vname, fi)
            fi.add_field((ptype, cname), function=NullFunc,
                         particle_type = True)
            df.append((ptype, cname))
            # Now we add some translations.
            df += self._setup_particle_type(ptype)
        self._derived_fields_add(df)

    def _setup_unknown_fields(self, list_of_fields = None, field_info = None,
                              skip_removal = False):
        field_info = field_info or self._fieldinfo_known
        field_list = list_of_fields or self.field_list
        dftype = self.default_fluid_type
        mylog.debug("Checking %s", field_list)
        for field in field_list:
            # By allowing a backup, we don't mandate that it's found in our
            # current field info.  This means we'll instead simply override
            # it.
            if not skip_removal:
                ff = self.field_info.pop(field, None)
            if field not in field_info:
                # Now we check if it's a gas field or what ...
                if isinstance(field, tuple) and field[0] in self.particle_types:
                    particle_type = True
                else:
                    particle_type = False
                if isinstance(field, tuple) and not particle_type and \
                   field[1] in field_info:
                    mylog.debug("Adding known field %s to list of fields", field)
                    self.field_info[field] = field_info[field[1]]
                    continue
                rootloginfo("Adding unknown field %s to list of fields", field)
                cf = None
                if self.has_key(field):
                    def external_wrapper(f):
                        def _convert_function(data):
                            return data.convert(f)
                        return _convert_function
                    cf = external_wrapper(field)
                # Note that we call add_field on the field_info directly.  This
                # will allow the same field detection mechanism to work for 1D, 2D
                # and 3D fields.
                self.field_info.add_field(
                        field, NullFunc, particle_type = particle_type,
                        convert_function=cf, take_log=False, units=r"Unknown")
            else:
                mylog.debug("Adding known field %s to list of fields", field)
                self.field_info[field] = field_info[field]

    def _setup_derived_fields(self):
        self.derived_field_list = []
        self.filtered_particle_types = []
        fc = self._derived_fields_to_check()
        self._derived_fields_add(fc)

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
            self._setup_particle_types([filter.name])
        return available

    def _derived_fields_to_check(self):
        fi = self.field_info
        # First we construct our list of fields to check
        fields_to_check = []
        for field in fi:
            finfo = fi[field]
            # Explicitly defined
            if isinstance(field, tuple):
                fields_to_check.append(field)
                continue
            # This one is implicity defined for all particle or fluid types.
            # So we check each.
            if not finfo.particle_type:
                fields_to_check.append(field)
                continue
            # We do a special case for 'all' later
            new_fields = []
            for pt in self.particle_types:
                new_fi = copy.copy(finfo)
                new_fi.name = (pt, new_fi.name)
                fi[new_fi.name] = new_fi
                new_fields.append(new_fi.name)
            fields_to_check += new_fields
        return fields_to_check

    def _derived_fields_add(self, fields_to_check = None):
        if fields_to_check is None:
            fields_to_check = []
        fi = self.field_info
        self._check_later = []
        for field in fields_to_check:
            try:
                fd = fi[field].get_dependencies(pf = self)
            except Exception as e:
                self._check_later.append(field)
                if type(e) != YTFieldNotFound:
                    mylog.debug("Raises %s during field %s detection.",
                                str(type(e)), field)
                continue
            missing = False
            # This next bit checks that we can't somehow generate everything.
            # We also manually update the 'requested' attribute
            requested = []
            for f in fd.requested:
                if (field[0], f) in self.field_list:
                    requested.append( (field[0], f) )
                elif f in self.field_list:
                    requested.append( f )
                elif isinstance(f, tuple) and f[1] in self.field_list:
                    requested.append( f )
                else:
                    missing = True
                    break
            if not missing: self.derived_field_list.append(field)
            fd.requested = set(requested)
            self.field_dependencies[field] = fd
            if not fi[field].particle_type and not isinstance(field, tuple):
                # Manually hardcode to 'gas'
                self.field_dependencies["gas", field] = fd
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    # Now all the object related stuff
    def all_data(self, find_max=False):
        if find_max: c = self.find_max("Density")[1]
        else: c = (self.domain_right_edge + self.domain_left_edge)/2.0
        return self.region(c,
            self.domain_left_edge, self.domain_right_edge)

    def _setup_data_io(self):
        if getattr(self, "io", None) is not None: return
        self.io = io_registry[self.dataset_type](self)

    def _split_fields(self, fields):
        # This will split fields into either generated or read fields
        fields_to_read, fields_to_generate = [], []
        for ftype, fname in fields:
            if fname in self.field_list or (ftype, fname) in self.field_list:
                fields_to_read.append((ftype, fname))
            else:
                fields_to_generate.append((ftype, fname))
        return fields_to_read, fields_to_generate

    def _read_particle_fields(self, fields, dobj, chunk = None):
        if len(fields) == 0: return {}, []
        selector = dobj.selector
        if chunk is None:
            self.index._identify_base_chunk(dobj)
        fields_to_return = {}
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_particle_selection(
            self.index._chunk_io(dobj, cache = False),
            selector,
            fields_to_read)
        for field in fields_to_read:
            ftype, fname = field
            finfo = self._get_field_info(*field)
            finfo = self._get_field_info(*field)
            conv_factor = finfo._convert_function(self)
            np.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        return fields_to_return, fields_to_generate

    def _read_fluid_fields(self, fields, dobj, chunk = None):
        if len(fields) == 0: return {}, []
        selector = dobj.selector
        if chunk is None:
            self.index._identify_base_chunk(dobj)
            chunk_size = dobj.size
        else:
            chunk_size = chunk.data_size
        fields_to_return = {}
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_fluid_selection(
            self.index._chunk_io(dobj, cache = False),
            selector,
            fields_to_read,
            chunk_size)
        for field in fields_to_read:
            ftype, fname = field
            conv_factor = self.field_info[fname]._convert_function(self)
            np.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        #mylog.debug("Don't know how to read %s", fields_to_generate)
        return fields_to_return, fields_to_generate

    def convert(self, unit):
        return self.conversion_factors[unit]

def _reconstruct_pf(*args, **kwargs):
    pfs = ParameterFileStore()
    pf = pfs.get_pf_hash(*args)
    return pf


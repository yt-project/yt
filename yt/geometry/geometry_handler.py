"""
Geometry container base class.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import cPickle
import weakref
import h5py
from exceptions import IOError, TypeError
from types import ClassType
import numpy as np
import abc
import copy

from yt.funcs import *
from yt.config import ytcfg
from yt.data_objects.data_containers import \
    data_object_registry
from yt.data_objects.field_info_container import \
    NullFunc
from yt.data_objects.particle_fields import \
    particle_deposition_functions
from yt.utilities.io_handler import io_registry
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_splitter
from yt.utilities.exceptions import YTFieldNotFound

class GeometryHandler(ParallelAnalysisInterface):
    _global_mesh = True

    def __init__(self, pf, data_style):
        ParallelAnalysisInterface.__init__(self)
        self.parameter_file = weakref.proxy(pf)
        self.pf = self.parameter_file

        self._initialize_state_variables()

        mylog.debug("Initializing data storage.")
        self._initialize_data_storage()

        # Must be defined in subclass
        mylog.debug("Setting up classes.")
        self._setup_classes()

        mylog.debug("Setting up domain geometry.")
        self._setup_geometry()

        mylog.debug("Initializing data grid data IO")
        self._setup_data_io()

        mylog.debug("Detecting fields.")
        self._detect_fields()

        mylog.debug("Detecting fields in backup.")
        self._detect_fields_backup()

        mylog.debug("Adding unknown detected fields")
        self._setup_unknown_fields()

        mylog.debug("Setting up derived fields")
        self._setup_derived_fields()

    def __del__(self):
        if self._data_file is not None:
            self._data_file.close()

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

    def _initialize_state_variables(self):
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}
        self.num_grids = None

    def _setup_classes(self, dd):
        # Called by subclass
        self.object_types = []
        self.objects = []
        self.plots = []
        for name, cls in sorted(data_object_registry.items()):
            cname = cls.__name__
            if cname.endswith("Base"): cname = cname[:-4]
            self._add_object_class(name, cname, cls, dd)
        if self.pf.refine_by != 2 and hasattr(self, 'proj') and \
            hasattr(self, 'overlap_proj'):
            mylog.warning("Refine by something other than two: reverting to"
                        + " overlap_proj")
            self.proj = self.overlap_proj
        if self.pf.dimensionality < 3 and hasattr(self, 'proj') and \
            hasattr(self, 'overlap_proj'):
            mylog.warning("Dimensionality less than 3: reverting to"
                        + " overlap_proj")
            self.proj = self.overlap_proj
        self.object_types.sort()

    def _setup_unknown_fields(self):
        known_fields = self.parameter_file._fieldinfo_known
        mylog.debug("Checking %s", self.field_list)
        for field in self.field_list:
            # By allowing a backup, we don't mandate that it's found in our
            # current field info.  This means we'll instead simply override
            # it.
            ff = self.parameter_file.field_info.pop(field, None)
            if field not in known_fields:
                if isinstance(field, types.TupleType) and \
                   field[0] in self.parameter_file.particle_types:
                    particle_type = True
                else:
                    particle_type = False
                rootloginfo("Adding unknown field %s to list of fields", field)
                cf = None
                if self.parameter_file.has_key(field):
                    def external_wrapper(f):
                        def _convert_function(data):
                            return data.convert(f)
                        return _convert_function
                    cf = external_wrapper(field)
                # Note that we call add_field on the field_info directly.  This
                # will allow the same field detection mechanism to work for 1D, 2D
                # and 3D fields.
                self.pf.field_info.add_field(
                        field, NullFunc, particle_type = particle_type,
                        convert_function=cf, take_log=False, units=r"Unknown")
            else:
                mylog.debug("Adding known field %s to list of fields", field)
                self.parameter_file.field_info[field] = known_fields[field]

    def _setup_derived_fields(self):
        self.derived_field_list = []
        self.filtered_particle_types = []
        fc, fac = self._derived_fields_to_check()
        self._derived_fields_add(fc, fac)

    def _setup_filtered_type(self, filter):
        if not filter.available(self.derived_field_list):
            return False
        fi = self.parameter_file.field_info
        fd = self.parameter_file.field_dependencies
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
            self.parameter_file.particle_types += (filter.name,)
            self.filtered_particle_types.append(filter.name)
            self._setup_particle_fields(filter.name, True)
        return available

    def _setup_particle_fields(self, ptype, filtered = False):
        pf = self.parameter_file
        pmass = self.parameter_file._particle_mass_name
        pcoord = self.parameter_file._particle_coordinates_name
        if pmass is None or pcoord is None: return
        df = particle_deposition_functions(ptype,
            pcoord, pmass, self.parameter_file.field_info)
        self._derived_fields_add(df)

    def _derived_fields_to_check(self):
        fi = self.parameter_file.field_info
        # First we construct our list of fields to check
        fields_to_check = []
        fields_to_allcheck = []
        for field in fi.keys():
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
            for pt in self.parameter_file.particle_types:
                new_fi = copy.copy(finfo)
                new_fi.name = (pt, new_fi.name)
                fi[new_fi.name] = new_fi
                new_fields.append(new_fi.name)
            fields_to_check += new_fields
            fields_to_allcheck.append(field)
        return fields_to_check, fields_to_allcheck

    def _derived_fields_add(self, fields_to_check = None,
                            fields_to_allcheck = None):
        if fields_to_check is None:
            fields_to_check = []
        if fields_to_allcheck is None:
            fields_to_allcheck = []
        fi = self.parameter_file.field_info
        for field in fields_to_check:
            try:
                fd = fi[field].get_dependencies(pf = self.parameter_file)
            except Exception as e:
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
            self.parameter_file.field_dependencies[field] = fd
            if not fi[field].particle_type and not isinstance(field, tuple):
                # Manually hardcode to 'gas'
                self.parameter_file.field_dependencies["gas", field] = fd
        for base_field in fields_to_allcheck:
            # Now we expand our field_info with the new fields
            all_available = all(((pt, field) in self.derived_field_list
                                 for pt in self.parameter_file.particle_types))
            if all_available:
                self.derived_field_list.append( ("all", field) )
                fi["all", base_field] = fi[base_field]
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    # Now all the object related stuff
    def all_data(self, find_max=False):
        pf = self.parameter_file
        if find_max: c = self.find_max("Density")[1]
        else: c = (pf.domain_right_edge + pf.domain_left_edge)/2.0
        return self.region(c,
            pf.domain_left_edge, pf.domain_right_edge)

    def _initialize_data_storage(self):
        if not ytcfg.getboolean('yt','serialize'): return
        fn = self.pf.storage_filename
        if fn is None:
            if os.path.isfile(os.path.join(self.directory,
                                "%s.yt" % self.pf.unique_identifier)):
                fn = os.path.join(self.directory,"%s.yt" % self.pf.unique_identifier)
            else:
                fn = os.path.join(self.directory,
                        "%s.yt" % self.parameter_file.basename)
        dir_to_check = os.path.dirname(fn)
        if dir_to_check == '':
            dir_to_check = '.'
        # We have four options:
        #    Writeable, does not exist      : create, open as append
        #    Writeable, does exist          : open as append
        #    Not writeable, does not exist  : do not attempt to open
        #    Not writeable, does exist      : open as read-only
        exists = os.path.isfile(fn)
        if not exists:
            writeable = os.access(dir_to_check, os.W_OK)
        else:
            writeable = os.access(fn, os.W_OK)
        writeable = writeable and not ytcfg.getboolean('yt','onlydeserialize')
        # We now have our conditional stuff
        self.comm.barrier()
        if not writeable and not exists: return
        if writeable:
            try:
                if not exists: self.__create_data_file(fn)
                self._data_mode = 'a'
            except IOError:
                self._data_mode = None
                return
        else:
            self._data_mode = 'r'

        self.__data_filename = fn
        self._data_file = h5py.File(fn, self._data_mode)

    def __create_data_file(self, fn):
        # Note that this used to be parallel_root_only; it no longer is,
        # because we have better logic to decide who owns the file.
        f = h5py.File(fn, 'a')
        f.close()

    def _setup_data_io(self):
        if getattr(self, "io", None) is not None: return
        self.io = io_registry[self.data_style]()

    def _save_data(self, array, node, name, set_attr=None, force=False, passthrough = False):
        """
        Arbitrary numpy data will be saved to the region in the datafile
        described by *node* and *name*.  If data file does not exist, it throws
        no error and simply does not save.
        """

        if self._data_mode != 'a': return
        try:
            node_loc = self._data_file[node]
            if name in node_loc and force:
                mylog.info("Overwriting node %s/%s", node, name)
                del self._data_file[node][name]
            elif name in node_loc and passthrough:
                return
        except:
            pass
        myGroup = self._data_file['/']
        for q in node.split('/'):
            if q: myGroup = myGroup.require_group(q)
        arr = myGroup.create_dataset(name,data=array)
        if set_attr is not None:
            for i, j in set_attr.items(): arr.attrs[i] = j
        self._data_file.flush()

    def _reload_data_file(self, *args, **kwargs):
        if self._data_file is None: return
        self._data_file.close()
        del self._data_file
        self._data_file = h5py.File(self.__data_filename, self._data_mode)

    save_data = parallel_splitter(_save_data, _reload_data_file)

    def _get_data_reader_dict(self):
        dd = { 'pf' : self.parameter_file, # Already weak
               'hierarchy': weakref.proxy(self) }
        return dd

    def _reset_save_data(self,round_robin=False):
        if round_robin:
            self.save_data = self._save_data
        else:
            self.save_data = parallel_splitter(self._save_data, self._reload_data_file)

    def save_object(self, obj, name):
        """
        Save an object (*obj*) to the data_file using the Pickle protocol,
        under the name *name* on the node /Objects.
        """
        s = cPickle.dumps(obj, protocol=-1)
        self.save_data(np.array(s, dtype='c'), "/Objects", name, force = True)

    def load_object(self, name):
        """
        Load and return and object from the data_file using the Pickle protocol,
        under the name *name* on the node /Objects.
        """
        obj = self.get_data("/Objects", name)
        if obj is None:
            return
        obj = cPickle.loads(obj.value)
        if iterable(obj) and len(obj) == 2:
            obj = obj[1] # Just the object, not the pf
        if hasattr(obj, '_fix_pickle'): obj._fix_pickle()
        return obj

    def get_data(self, node, name):
        """
        Return the dataset with a given *name* located at *node* in the
        datafile.
        """
        if self._data_file == None:
            return None
        if node[0] != "/": node = "/%s" % node

        myGroup = self._data_file['/']
        for group in node.split('/'):
            if group:
                if group not in myGroup:
                    return None
                myGroup = myGroup[group]
        if name not in myGroup:
            return None

        full_name = "%s/%s" % (node, name)
        try:
            return self._data_file[full_name][:]
        except TypeError:
            return self._data_file[full_name]

    def _close_data_file(self):
        if self._data_file:
            self._data_file.close()
            del self._data_file
            self._data_file = None

    def _add_object_class(self, name, class_name, base, dd):
        self.object_types.append(name)
        obj = type(class_name, (base,), dd)
        setattr(self, name, obj)

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
            self._identify_base_chunk(dobj)
        fields_to_return = {}
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_particle_selection(
            self._chunk_io(dobj, cache = False),
            selector,
            fields_to_read)
        for field in fields_to_read:
            ftype, fname = field
            finfo = self.pf._get_field_info(*field)
            conv_factor = finfo._convert_function(self)
            np.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        return fields_to_return, fields_to_generate

    def _read_fluid_fields(self, fields, dobj, chunk = None):
        if len(fields) == 0: return {}, []
        selector = dobj.selector
        if chunk is None:
            self._identify_base_chunk(dobj)
            chunk_size = dobj.size
        else:
            chunk_size = chunk.data_size
        fields_to_return = {}
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_fluid_selection(
            self._chunk_io(dobj, cache = False),
            selector,
            fields_to_read,
            chunk_size)
        for field in fields_to_read:
            ftype, fname = field
            conv_factor = self.pf.field_info[fname]._convert_function(self)
            np.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        #mylog.debug("Don't know how to read %s", fields_to_generate)
        return fields_to_return, fields_to_generate


    def _chunk(self, dobj, chunking_style, ngz = 0, **kwargs):
        # A chunk is either None or (grids, size)
        if dobj._current_chunk is None:
            self._identify_base_chunk(dobj)
        if ngz != 0 and chunking_style != "spatial":
            raise NotImplementedError
        if chunking_style == "all":
            return self._chunk_all(dobj, **kwargs)
        elif chunking_style == "spatial":
            return self._chunk_spatial(dobj, ngz, **kwargs)
        elif chunking_style == "io":
            return self._chunk_io(dobj, **kwargs)
        else:
            raise NotImplementedError

def cached_property(func):
    n = '_%s' % func.func_name
    def cached_func(self):
        if self._cache and getattr(self, n, None) is not None:
            return getattr(self, n)
        if self.data_size is None:
            tr = self._accumulate_values(n[1:])
        else:
            tr = func(self)
        if self._cache:
            setattr(self, n, tr)
        return tr
    return property(cached_func)

class YTDataChunk(object):

    def __init__(self, dobj, chunk_type, objs, data_size = None,
                 field_type = None, cache = False):
        self.dobj = dobj
        self.chunk_type = chunk_type
        self.objs = objs
        self.data_size = data_size
        self._field_type = field_type
        self._cache = cache

    def _accumulate_values(self, method):
        # We call this generically.  It's somewhat slower, since we're doing
        # costly getattr functions, but this allows us to generalize.
        mname = "select_%s" % method
        arrs = []
        for obj in self.objs:
            f = getattr(obj, mname)
            arrs.append(f(self.dobj))
        arrs = np.concatenate(arrs)
        self.data_size = arrs.shape[0]
        return arrs

    @cached_property
    def fcoords(self):
        ci = np.empty((self.data_size, 3), dtype='float64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_fcoords(self.dobj)
            if c.shape[0] == 0: continue
            ci[ind:ind+c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cached_property
    def icoords(self):
        ci = np.empty((self.data_size, 3), dtype='int64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_icoords(self.dobj)
            if c.shape[0] == 0: continue
            ci[ind:ind+c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cached_property
    def fwidth(self):
        ci = np.empty((self.data_size, 3), dtype='float64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_fwidth(self.dobj)
            if c.shape[0] == 0: continue
            ci[ind:ind+c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cached_property
    def ires(self):
        ci = np.empty(self.data_size, dtype='int64')
        if self.data_size == 0: return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_ires(self.dobj)
            if c.shape == 0: continue
            ci[ind:ind+c.size] = c
            ind += c.size
        return ci

    @cached_property
    def tcoords(self):
        self.dtcoords
        return self._tcoords

    @cached_property
    def dtcoords(self):
        ct = np.empty(self.data_size, dtype='float64')
        cdt = np.empty(self.data_size, dtype='float64')
        self._tcoords = ct # Se this for tcoords
        if self.data_size == 0: return cdt
        ind = 0
        for obj in self.objs:
            gdt, gt = obj.tcoords(self.dobj)
            if gt.shape == 0: continue
            ct[ind:ind+gt.size] = gt
            cdt[ind:ind+gdt.size] = gdt
            ind += gt.size
        return cdt

class ChunkDataCache(object):
    def __init__(self, base_iter, preload_fields, geometry_handler,
                 max_length = 256):
        # At some point, max_length should instead become a heuristic function,
        # potentially looking at estimated memory usage.  Note that this never
        # initializes the iterator; it assumes the iterator is already created,
        # and it calls next() on it.
        self.base_iter = base_iter.__iter__()
        self.queue = []
        self.max_length = max_length
        self.preload_fields = preload_fields
        self.geometry_handler = geometry_handler
        self.cache = {}

    def __iter__(self):
        return self
    
    def next(self):
        if len(self.queue) == 0:
            for i in range(self.max_length):
                try:
                    self.queue.append(self.base_iter.next())
                except StopIteration:
                    break
            # If it's still zero ...
            if len(self.queue) == 0: raise StopIteration
            chunk = YTDataChunk(None, "cache", self.queue, cache=False)
            self.cache = self.geometry_handler.io._read_chunk_data(
                chunk, self.preload_fields)
        g = self.queue.pop(0)
        g._initialize_cache(self.cache.pop(g.id, {}))
        return g

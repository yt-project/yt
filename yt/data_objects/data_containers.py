"""
The base classes for selecting and returning data.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Author: Britton Smith <Britton.Smith@colorado.edu>
Affiliation: University of Colorado at Boulder
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

import itertools
import types

data_object_registry = {}

import numpy as na
import weakref
import shelve
from contextlib import contextmanager

from yt.funcs import *

from yt.data_objects.particle_io import particle_handler_registry
from yt.utilities.amr_utils import \
    march_cubes_grid, march_cubes_grid_flux
from yt.utilities.definitions import  x_dict, y_dict
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.parameter_file_storage import \
    ParameterFileStore
from .derived_quantities import DerivedQuantityCollection
from .field_info_container import \
    NeedsGridType
import yt.geometry.selection_routines

def force_array(item, shape):
    try:
        sh = item.shape
        return item
    except AttributeError:
        if item:
            return na.ones(shape, dtype='bool')
        else:
            return na.zeros(shape, dtype='bool')

class YTFieldData(dict):
    """
    A Container object for field data, instead of just having it be a dict.
    """
    pass
        

class YTDataContainer(object):
    """
    Generic YTDataContainer container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    and deal with passing back and forth field parameters.
    """
    _grids = _domains = None
    _num_ghost_zones = 0
    _con_args = ()
    _skip_add = False
    _container_fields = ()

    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_type_name") and not cls._skip_add:
                data_object_registry[cls._type_name] = cls

    def __init__(self, pf, fields, **kwargs):
        """
        Typically this is never called directly, but only due to inheritance.
        It associates a :class:`~yt.data_objects.api.StaticOutput` with the class,
        sets its initial set of fields, and the remainder of the arguments
        are passed as field_parameters.
        """
        if pf != None:
            self.pf = pf
            self.hierarchy = pf.hierarchy
        self.hierarchy.objects.append(weakref.proxy(self))
        mylog.debug("Appending object to %s (type: %s)", self.pf, type(self))
        if fields == None: fields = []
        self.fields = ensure_list(fields)[:]
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.__set_default_field_parameters()
        self._cut_masks = {}
        self._point_indices = {}
        self._vc_data = {}
        for key, val in kwargs.items():
            mylog.debug("Setting %s to %s", key, val)
            self.set_field_parameter(key, val)

    def __set_default_field_parameters(self):
        self.set_field_parameter("center",na.zeros(3,dtype='float64'))
        self.set_field_parameter("bulk_velocity",na.zeros(3,dtype='float64'))

    def _set_center(self, center):
        if center is None:
            pass
        elif isinstance(center, (types.ListType, types.TupleType, na.ndarray)):
            center = na.array(center)
        elif center in ("center", "c"): # is this dangerous for race conditions?
            center = self.pf.domain_center
        elif center == "max": # is this dangerous for race conditions?
            center = self.pf.h.find_max("Density")[1]
        elif center.startswith("max_"):
            center = self.pf.h.find_max(center[4:])[1]
        else:
            center = na.array(center, dtype='float64')
        self.center = center
        self.set_field_parameter('center', center)

    def get_field_parameter(self, name, default=None):
        """
        This is typically only used by derived field functions, but
        it returns parameters used to generate fields.
        """
        if self.field_parameters.has_key(name):
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
        return self.field_parameters.has_key(name)

    def convert(self, datatype):
        """
        This will attempt to convert a given unit to cgs from code units.
        It either returns the multiplicative factor or throws a KeyError.
        """
        return self.pf[datatype]

    def clear_data(self):
        """
        Clears out all data from the YTDataContainer instance, freeing memory.
        """
        self.field_data.clear()

    def clear_cache(self):
        """
        Clears out all cache, freeing memory.
        """
        for _cm in self._cut_masks: del _cm
        for _pi in self._point_indices: del _pi
        for _field in self._vc_data:
            for _vc in _field: del _vc

    def has_key(self, key):
        """
        Checks if a data field already exists.
        """
        return self.field_data.has_key(key)

    def keys(self):
        return self.field_data.keys()

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        if key not in self.field_data:
            self.get_data(key)
        f = self._determine_fields(key)[0]
        return self.field_data[f]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        if key not in self.fields: self.fields.append(key)
        self.field_data[key] = val

    def __delitem__(self, key):
        """
        Deletes a field
        """
        try:
            del self.fields[self.fields.index(key)]
        except ValueError:
            pass
        del self.field_data[key]

    def _generate_field(self, field):
        ftype, fname = field
        if fname in self._container_fields:
            return self._generate_container_field(field)
        elif fname not in self.pf.field_info:
            raise KeyError(field)
        elif self.pf.field_info[fname].particle_type:
            return self._generate_particle_field(field)
        else:
            return self._generate_fluid_field(field)

    def _generate_fluid_field(self, field):
        # First we check the validator
        ftype, fname = field
        if self._current_chunk is None or \
           self._current_chunk.chunk_type != "spatial":
            gen_obj = self
        else:
            gen_obj = self._current_chunk.objs[0]
        try:
            self.pf.field_info[fname].check_available(gen_obj)
        except NeedsGridType, ngt_exception:
            rv = na.empty(self.size, dtype="float64")
            ind = 0
            ngz = ngt_exception.ghost_zones
            for io_chunk in self.chunks([], "io"):
                for i,chunk in enumerate(self.chunks(field, "spatial", ngz = ngz)):
                    mask = self._current_chunk.objs[0].select(self.selector)
                    if mask is None: continue
                    data = self[field]
                    if ngz > 0:
                        data = data[ngz:-ngz, ngz:-ngz, ngz:-ngz]
                    data = data[mask]
                    rv[ind:ind+data.size] = data
                    ind += data.size
        else:
            rv = self.pf.field_info[fname](gen_obj)
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
            self.pf.field_info[fname].check_available(gen_obj)
        except NeedsGridType, ngt_exception:
            if ngt_exception.ghost_zones != 0:
                raise NotImplementedError
            size = self._count_particles(ftype)
            rv = na.empty(size, dtype="float64")
            ind = 0
            for io_chunk in self.chunks([], "io"):
                for i,chunk in enumerate(self.chunks(field, "spatial")):
                    x, y, z = (self[ftype, 'particle_position_%s' % ax]
                               for ax in 'xyz')
                    mask = self._current_chunk.objs[0].select_particles(
                        self.selector, x, y, z)
                    if mask is None: continue
                    # This requests it from the grid and does NOT mask it
                    data = self[field][mask]
                    rv[ind:ind+data.size] = data
                    ind += data.size
        else:
            rv = self.pf.field_info[fname](gen_obj)
        return rv

    def _count_particles(self, ftype):
        for (f1, f2), val in self.field_data.items():
            if f1 == ftype:
                return val.size
        size = 0
        for io_chunk in self.chunks([], "io"):
            for i,chunk in enumerate(self.chunks([], "spatial")):
                x, y, z = (self[ftype, 'particle_position_%s' % ax]
                            for ax in 'xyz')
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
        field_data = na.array([self.field_data[field] for field in field_order])
        for line in range(field_data.shape[1]):
            field_data[:,line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

    def save_object(self, name, filename = None):
        """
        Save an object.  If *filename* is supplied, it will be stored in
        a :mod:`shelve` file of that name.  Otherwise, it will be stored via
        :meth:`yt.data_objects.api.GridGeometryHandler.save_object`.
        """
        if filename is not None:
            ds = shelve.open(filename, protocol=-1)
            if name in ds:
                mylog.info("Overwriting %s in %s", name, filename)
            ds[name] = self
            ds.close()
        else:
            self.hierarchy.save_object(self, name)

    def __reduce__(self):
        args = tuple([self.pf._hash(), self._type_name] +
                     [getattr(self, n) for n in self._con_args] +
                     [self.field_parameters])
        return (_reconstruct_object, args)

    def __repr__(self, clean = False):
        # We'll do this the slow way to be clear what's going on
        if clean: s = "%s: " % (self.__class__.__name__)
        else: s = "%s (%s): " % (self.__class__.__name__, self.pf)
        s += ", ".join(["%s=%s" % (i, getattr(self,i))
                       for i in self._con_args])
        return s

    def _get_field_info(self, fname):
        if fname not in self.pf.field_info:
            raise YTFieldNotFound(fname, self.pf)
        finfo = self.pf.field_info[fname]
        return finfo

    def _determine_fields(self, fields):
        fields = ensure_list(fields)
        explicit_fields = []
        for field in fields:
            if isinstance(field, types.TupleType):
                if len(field) != 2 or \
                   not isinstance(field[0], types.StringTypes) or \
                   not isinstance(field[1], types.StringTypes):
                    raise YTFieldNotParseable(field)
                ftype, fname = field
                finfo = self._get_field_info(fname)
                explicit_fields.append(field)
            else:
                fname = field
                finfo = self._get_field_info(fname)
                if finfo.particle_type:
                    ftype = "all"
                else:
                    ftype = self.pf.default_fluid_type
            if finfo.particle_type and ftype not in self.pf.particle_types:
                raise YTFieldTypeNotFound(ftype)
            elif not finfo.particle_type and ftype not in self.pf.fluid_types:
                raise YTFieldTypeNotFound(ftype)
            explicit_fields.append((ftype, fname))
        return explicit_fields

class GenerationInProgress(Exception):
    def __init__(self, fields):
        self.fields = fields
        super(GenerationInProgress, self).__init__()

class YTSelectionContainer(YTDataContainer, ParallelAnalysisInterface):
    _locked = False
    _sort_by = None
    _selector = None
    _current_chunk = None
    size = None
    shape = None

    def __init__(self, *args, **kwargs):
        super(YTSelectionContainer, self).__init__(*args, **kwargs)

    @property
    def selector(self):
        if self._selector is not None: return self._selector
        sclass = getattr(yt.geometry.selection_routines,
                         "%s_selector" % self._type_name, None)
        if sclass is None: raise NotImplementedError
        self._selector = sclass(self)
        return self._selector

    def chunks(self, fields, chunking_style, **kwargs):
        # This is an iterator that will yield the necessary chunks.
        self.get_data() # Ensure we have built ourselves
        if fields is None: fields = []
        for chunk in self.hierarchy._chunk(self, chunking_style, **kwargs):
            with self._chunked_read(chunk):
                self.get_data(fields)
                # NOTE: we yield before releasing the context
                yield self

    def get_data(self, fields=None):
        if self._current_chunk is None:
            self.hierarchy._identify_base_chunk(self)
        if fields is None: return
        fields = self._determine_fields(fields)
        # Now we collect all our fields
        fields_to_get = [f for f in fields if f not in self.field_data]
        if len(fields_to_get) == 0:
            return
        elif self._locked == True:
            raise GenerationInProgress(fields)
        # At this point, we want to figure out *all* our dependencies.
        inspected = 0
        for ftype, field in itertools.cycle(fields_to_get):
            if inspected >= len(fields_to_get): break
            inspected += 1
            if field not in self.pf.field_dependencies: continue
            fd = self.pf.field_dependencies[field]
            requested = self._determine_fields(fd.requested)
            deps = [d for d in requested if d not in fields_to_get]
            fields_to_get += deps
        # We now split up into readers for the types of fields
        fluids, particles = [], []
        for ftype, fname in fields_to_get:
            if self.pf.field_info[fname].particle_type:
                particles.append((ftype, fname))
            else:
                fluids.append((ftype, fname))
        # The _read method will figure out which fields it needs to get from
        # disk, and return a dict of those fields along with the fields that
        # need to be generated.
        read_fluids, gen_fluids = self.hierarchy._read_fluid_fields(
                                        fluids, self, self._current_chunk)
        self.field_data.update(read_fluids)

        read_particles, gen_particles = self.hierarchy._read_particle_fields(
                                        particles, self, self._current_chunk)
        self.field_data.update(read_particles)
        fields_to_generate = gen_fluids + gen_particles
        index = 0
        with self._field_lock():
            while any(f not in self.field_data for f in fields_to_generate):
                field = fields_to_generate[index % len(fields_to_generate)]
                index += 1
                if field in self.field_data: continue
                try:
                    self.field_data[field] = self._generate_field(field)
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
        old_size, self.size = self.size, chunk.data_size
        old_chunk, self._current_chunk = self._current_chunk, chunk
        old_locked, self._locked = self._locked, False
        self.shape = (self.size,)
        yield
        self.field_data = old_field_data
        self.size = old_size
        self.shape = (old_size,)
        self._current_chunk = old_chunk
        self._locked = old_locked

    @property
    def icoords(self):
        return self._current_chunk.icoords

    @property
    def fcoords(self):
        return self._current_chunk.fcoords

    @property
    def ires(self):
        return self._current_chunk.ires

    @property
    def fwidth(self):
        return self._current_chunk.fwidth

class YTSelectionContainer1D(YTSelectionContainer):
    _spatial = False
    def __init__(self, pf, fields, **kwargs):
        super(YTSelectionContainer1D, self).__init__(pf, fields, **kwargs)
        self._grids = None
        self._sortkey = None
        self._sorted = {}

class YTSelectionContainer2D(YTSelectionContainer):
    _key_fields = ['px','py','pdx','pdy']
    """
    Class to represent a set of :class:`YTDataContainer` that's 2-D in nature, and
    thus does not have as many actions as the 3-D data types.
    """
    _spatial = False
    def __init__(self, axis, fields, pf=None, **kwargs):
        """
        Prepares the YTSelectionContainer2D, normal to *axis*.  If *axis* is 4, we are not
        aligned with any axis.
        """
        ParallelAnalysisInterface.__init__(self)
        self.axis = axis
        super(YTSelectionContainer2D, self).__init__(pf, fields, **kwargs)
        self.field = ensure_list(fields)[0]
        self.set_field_parameter("axis",axis)
        
    def _convert_field_name(self, field):
        return field

    def to_frb(self, width, resolution, center = None):
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
        resolution : int or tuple of ints
            The number of pixels on a side of the final FRB.
        center : array-like of floats, optional
            The center of the FRB.  If not specified, defaults to the center of
            the current object.

        Returns
        -------
        frb : :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`
            A fixed resolution buffer, which can be queried for fields.

        Examples
        --------

        >>> proj = pf.h.proj(0, "Density")
        >>> frb = proj.to_frb( (100.0, 'kpc'), 1024)
        >>> write_image(na.log10(frb["Density"]), 'density_100kpc.png')
        """
        if center is None:
            center = self.get_field_parameter("center")
            if center is None:
                center = (self.pf.domain_right_edge
                        + self.pf.domain_left_edge)/2.0
        if iterable(width):
            w, u = width
            width = w/self.pf[u]
        if not iterable(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import FixedResolutionBuffer
        xax = x_dict[self.axis]
        yax = y_dict[self.axis]
        bounds = (center[xax] - width/2.0, center[xax] + width/2.0,
                  center[yax] - width/2.0, center[yax] + width/2.0)
        frb = FixedResolutionBuffer(self, bounds, resolution)
        return frb

    _okay_to_serialize = True

    def _store_fields(self, fields, node_name = None, force = False):
        fields = ensure_list(fields)
        if node_name is None: node_name = self._gen_node_name()
        for field in fields:
            #mylog.debug("Storing %s in node %s",
                #self._convert_field_name(field), node_name)
            self.hierarchy.save_data(self[field], node_name,
                self._convert_field_name(field), force = force,
                passthrough = True)

    def _obtain_fields(self, fields, node_name = None):
        if not self._okay_to_serialize: return
        fields = ensure_list(fields)
        if node_name is None: node_name = self._gen_node_name()
        for field in fields:
            #mylog.debug("Trying to obtain %s from node %s",
                #self._convert_field_name(field), node_name)
            fdata=self.hierarchy.get_data(node_name, 
                self._convert_field_name(field))
            if fdata is not None:
                #mylog.debug("Got %s from node %s", field, node_name)
                self[field] = fdata[:]
        return True

    def _deserialize(self, node_name = None):
        if not self._okay_to_serialize: return
        self._obtain_fields(self._key_fields, node_name)
        self._obtain_fields(self.fields, node_name)

    def _serialize(self, node_name = None, force = False):
        if not self._okay_to_serialize: return
        self._store_fields(self._key_fields, node_name, force)
        self._store_fields(self.fields, node_name, force)

class YTSelectionContainer3D(YTSelectionContainer):
    _key_fields = ['x','y','z','dx','dy','dz']
    """
    Class describing a cluster of data points, not necessarily sharing any
    particular attribute.
    """
    _spatial = False
    _num_ghost_zones = 0
    def __init__(self, center, fields, pf = None, **kwargs):
        """
        Returns an instance of YTSelectionContainer3D, or prepares one.  Usually only
        used as a base class.  Note that *center* is supplied, but only used
        for fields and quantities that require it.
        """
        ParallelAnalysisInterface.__init__(self)
        super(YTSelectionContainer3D, self).__init__(pf, fields, **kwargs)
        self._set_center(center)
        self.coords = None
        self._grids = None
        self.quantities = DerivedQuantityCollection(self)

    def cut_region(self, field_cuts):
        """
        Return an InLineExtractedRegion, where the grid cells are cut on the
        fly with a set of field_cuts.
        """
        return YTValueCutExtractionBase(self, field_cuts)

    def extract_region(self, indices):
        """
        Return an ExtractedRegion where the points contained in it are defined
        as the points in `this` data object with the given *indices*.
        """
        fp = self.field_parameters.copy()
        return YTSelectedIndicesBase(self, indices, **fp)

    def extract_isocontours(self, field, value, filename = None,
                            rescale = False, sample_values = None):
        r"""This identifies isocontours on a cell-by-cell basis, with no
        consideration of global connectedness, and returns the vertices of the
        Triangles in that isocontour.

        This function simply returns the vertices of all the triangles
        calculated by the marching cubes algorithm; for more complex
        operations, such as identifying connected sets of cells above a given
        threshold, see the extract_connected_sets function.  This is more
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

        References
        ----------

        .. [1] Marching Cubes: http://en.wikipedia.org/wiki/Marching_cubes

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        output the vertices to "triangles.obj" after rescaling them.

        >>> dd = pf.h.all_data()
        >>> rho = dd.quantities["WeightedAverageQuantity"](
        ...     "Density", weight="CellMassMsun")
        >>> verts = dd.extract_isocontours("Density", rho,
        ...             "triangles.obj", True)
        """
        verts = []
        samples = []
        for i, g in enumerate(self._get_grid_objs()):
            my_verts = self._extract_isocontours_from_grid(
                            g, field, value, sample_values)
            if sample_values is not None:
                my_verts, svals = my_verts
                samples.append(svals)
            verts.append(my_verts)
        verts = na.concatenate(verts).transpose()
        verts = self.comm.par_combine_object(verts, op='cat', datatype='array')
        verts = verts.transpose()
        if sample_values is not None:
            samples = na.concatenate(samples)
            samples = self.comm.par_combine_object(samples, op='cat',
                                datatype='array')
        if rescale:
            mi = na.min(verts, axis=0)
            ma = na.max(verts, axis=0)
            verts = (verts - mi) / (ma - mi).max()
        if filename is not None and self.comm.rank == 0:
            f = open(filename, "w")
            for v1 in verts:
                f.write("v %0.16e %0.16e %0.16e\n" % (v1[0], v1[1], v1[2]))
            for i in range(len(verts)/3):
                f.write("f %s %s %s\n" % (i*3+1, i*3+2, i*3+3))
            f.close()
        if sample_values is not None:
            return verts, samples
        return verts


    def _extract_isocontours_from_grid(self, grid, field, value,
                                       sample_values = None):
        mask = self._get_cut_mask(grid) * grid.child_mask
        vals = grid.get_vertex_centered_data(field, no_ghost = False)
        if sample_values is not None:
            svals = grid.get_vertex_centered_data(sample_values)
        else:
            svals = None
        my_verts = march_cubes_grid(value, vals, mask, grid.LeftEdge,
                                    grid.dds, svals)
        return my_verts

    def calculate_isocontour_flux(self, field, value,
                    field_x, field_y, field_z, fluxing_field = None):
        r"""This identifies isocontours on a cell-by-cell basis, with no
        consideration of global connectedness, and calculates the flux over
        those contours.

        This function will conduct marching cubes on all the cells in a given
        data container (grid-by-grid), and then for each identified triangular
        segment of an isocontour in a given cell, calculate the gradient (i.e.,
        normal) in the isocontoured field, interpolate the local value of the
        "fluxing" field, the area of the triangle, and then return:

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

        References
        ----------

        .. [1] Marching Cubes: http://en.wikipedia.org/wiki/Marching_cubes

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        calculate the metal flux over it.

        >>> dd = pf.h.all_data()
        >>> rho = dd.quantities["WeightedAverageQuantity"](
        ...     "Density", weight="CellMassMsun")
        >>> flux = dd.calculate_isocontour_flux("Density", rho,
        ...     "x-velocity", "y-velocity", "z-velocity", "Metal_Density")
        """
        flux = 0.0
        for g in self._get_grid_objs():
            flux += self._calculate_flux_in_grid(g, field, value,
                    field_x, field_y, field_z, fluxing_field)
        flux = self.comm.mpi_allreduce(flux, op="sum")
        return flux

    def _calculate_flux_in_grid(self, grid, field, value,
                    field_x, field_y, field_z, fluxing_field = None):
        mask = self._get_cut_mask(grid) * grid.child_mask
        vals = grid.get_vertex_centered_data(field)
        if fluxing_field is None:
            ff = na.ones(vals.shape, dtype="float64")
        else:
            ff = grid.get_vertex_centered_data(fluxing_field)
        xv, yv, zv = [grid.get_vertex_centered_data(f) for f in 
                     [field_x, field_y, field_z]]
        return march_cubes_grid_flux(value, vals, xv, yv, zv,
                    ff, mask, grid.LeftEdge, grid.dds)

    def extract_connected_sets(self, field, num_levels, min_val, max_val,
                                log_space=True, cumulative=True, cache=False):
        """
        This function will create a set of contour objects, defined
        by having connected cell structures, which can then be
        studied and used to 'paint' their source grids, thus enabling
        them to be plotted.
        """
        if log_space:
            cons = na.logspace(na.log10(min_val),na.log10(max_val),
                               num_levels+1)
        else:
            cons = na.linspace(min_val, max_val, num_levels+1)
        contours = {}
        if cache: cached_fields = defaultdict(lambda: dict())
        else: cached_fields = None
        for level in range(num_levels):
            contours[level] = {}
            if cumulative:
                mv = max_val
            else:
                mv = cons[level+1]
            from yt.analysis_modules.level_sets.api import identify_contours
            cids = identify_contours(self, field, cons[level], mv,
                                     cached_fields)
            for cid, cid_ind in cids.items():
                contours[level][cid] = self.extract_region(cid_ind)
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
            if default_value != None:
                grid[field] = na.ones(grid.ActiveDimensions)*default_value
            grid[field][self._get_point_indices(grid)] = value

    _particle_handler = None

    @property
    def particles(self):
        if self._particle_handler is None:
            self._particle_handler = \
                particle_handler_registry[self._type_name](self.pf, self)
        return self._particle_handler


    def volume(self, unit = "unitary"):
        """
        Return the volume of the data container in units *unit*.
        This is found by adding up the volume of the cells with centers
        in the container, rather than using the geometric shape of
        the container, so this may vary very slightly
        from what might be expected from the geometric volume.
        """
        return self.quantities["TotalQuantity"]("CellVolume")[0] * \
            (self.pf[unit] / self.pf['cm']) ** 3.0

def _reconstruct_object(*args, **kwargs):
    pfid = args[0]
    dtype = args[1]
    field_parameters = args[-1]
    # will be much nicer when we can do pfid, *a, fp = args
    args, new_args = args[2:-1], []
    for arg in args:
        if iterable(arg) and len(arg) == 2 \
           and not isinstance(arg, types.DictType) \
           and isinstance(arg[1], YTDataContainer):
            new_args.append(arg[1])
        else: new_args.append(arg)
    pfs = ParameterFileStore()
    pf = pfs.get_pf_hash(pfid)
    cls = getattr(pf.h, dtype)
    obj = cls(*new_args)
    obj.field_parameters.update(field_parameters)
    return pf, obj


class YTSelectedIndicesBase(YTSelectionContainer3D):
    """
    ExtractedRegions are arbitrarily defined containers of data, useful
    for things like selection along a baryon field.
    """
    _type_name = "extracted_region"
    _con_args = ('_base_region', '_indices')
    def __init__(self, base_region, indices, force_refresh=True, **kwargs):
        """An arbitrarily defined data container that allows for selection
        of all data meeting certain criteria.

        In order to create an arbitrarily selected set of data, the
        ExtractedRegion takes a `base_region` and a set of `indices`
        and creates a region within the `base_region` consisting of
        all data indexed by the `indices`. Note that `indices` must be
        precomputed. This does not work well for parallelized
        operations.

        Parameters
        ----------
        base_region : yt data source
            A previously selected data source.
        indices : array_like
            An array of indices

        Other Parameters
        ----------------
        force_refresh : bool
           Force a refresh of the data. Defaults to True.

        Examples
        --------
        """
        cen = kwargs.pop("center", None)
        if cen is None: cen = base_region.get_field_parameter("center")
        YTSelectionContainer3D.__init__(self, center=cen,
                            fields=None, pf=base_region.pf, **kwargs)
        self._base_region = base_region # We don't weakly reference because
                                        # It is not cyclic
        if isinstance(indices, types.DictType):
            self._indices = indices
            self._grids = self._base_region.pf.h.grids[self._indices.keys()]
        else:
            self._grids = None
            self._base_indices = indices

    def _get_cut_particle_mask(self, grid):
        # Override to provide a warning
        mylog.warning("Returning all particles from an Extracted Region.  This could be incorrect!")
        return True

    def _get_list_of_grids(self):
        # Okay, so what we're going to want to do is get the pointI from
        # region._get_point_indices(grid) for grid in base_region._grids,
        # and then construct an array of those, which we will select along indices.
        if self._grids != None: return
        grid_vals, xi, yi, zi = [], [], [], []
        for grid in self._base_region._grids:
            xit,yit,zit = self._base_region._get_point_indices(grid)
            grid_vals.append(na.ones(xit.shape, dtype='int') * (grid.id-grid._id_offset))
            xi.append(xit)
            yi.append(yit)
            zi.append(zit)
        grid_vals = na.concatenate(grid_vals)[self._base_indices]
        grid_order = na.argsort(grid_vals)
        # Note: grid_vals is still unordered
        grid_ids = na.unique(grid_vals)
        xi = na.concatenate(xi)[self._base_indices][grid_order]
        yi = na.concatenate(yi)[self._base_indices][grid_order]
        zi = na.concatenate(zi)[self._base_indices][grid_order]
        bc = na.bincount(grid_vals)
        splits = []
        for i,v in enumerate(bc):
            if v > 0: splits.append(v)
        splits = na.add.accumulate(splits)
        xis, yis, zis = [na.array_split(aa, splits) for aa in [xi,yi,zi]]
        self._indices = {}
        h = self._base_region.pf.h
        for grid_id, x, y, z in itertools.izip(grid_ids, xis, yis, zis):
            # grid_id needs no offset
            ll = h.grids[grid_id].ActiveDimensions.prod() \
               - (na.logical_not(h.grids[grid_id].child_mask)).sum()
            # This means we're completely enclosed, except for child masks
            if x.size == ll:
                self._indices[grid_id] = None
            else:
                # This will slow things down a bit, but conserve memory
                self._indices[grid_id] = \
                    na.zeros(h.grids[grid_id].ActiveDimensions, dtype='bool')
                self._indices[grid_id][(x,y,z)] = True
        self._grids = h.grids[self._indices.keys()]

    def _is_fully_enclosed(self, grid):
        if self._indices[grid.id-grid._id_offset] is None or \
            (self._indices[grid.id-grid._id_offset][0].size ==
             grid.ActiveDimensions.prod()):
            return True
        return False

    def _get_cut_mask(self, grid):
        cm = na.zeros(grid.ActiveDimensions, dtype='bool')
        cm[self._get_point_indices(grid, False)] = True
        return cm

    __empty_array = na.array([], dtype='bool')
    def _get_point_indices(self, grid, use_child_mask=True):
        # Yeah, if it's not true, we don't care.
        tr = self._indices.get(grid.id-grid._id_offset, self.__empty_array)
        if tr is None: tr = na.where(grid.child_mask)
        else: tr = na.where(tr)
        return tr

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = "%s (%s): " % (self.__class__.__name__, self.pf)
        s += ", ".join(["%s=%s" % (i, getattr(self,i))
                       for i in self._con_args if i != "_indices"])
        return s

    def join(self, other):
        ng = {}
        gs = set(self._indices.keys() + other._indices.keys())
        for g in gs:
            grid = self.pf.h.grids[g]
            if g in other._indices and g in self._indices:
                # We now join the indices
                ind = na.zeros(grid.ActiveDimensions, dtype='bool')
                ind[self._indices[g]] = True
                ind[other._indices[g]] = True
                if ind.prod() == grid.ActiveDimensions.prod(): ind = None
            elif g in self._indices:
                ind = self._indices[g]
            elif g in other._indices:
                ind = other._indices[g]
            # Okay we have indices
            if ind is not None: ind = ind.copy()
            ng[g] = ind
        gl = self.pf.h.grids[list(gs)]
        gc = self.pf.h.grid_collection(
            self._base_region.get_field_parameter("center"), gl)
        return self.pf.h.extracted_region(gc, ng)


class YTValueCutExtractionBase(YTSelectionContainer3D):
    """
    In-line extracted regions accept a base region and a set of field_cuts to
    determine which points in a grid should be included.
    """
    def __init__(self, base_region, field_cuts, **kwargs):
        cen = base_region.get_field_parameter("center")
        YTSelectionContainer3D.__init__(self, center=cen,
                            fields=None, pf=base_region.pf, **kwargs)
        self._base_region = base_region # We don't weakly reference because
                                        # It is not cyclic
        self._field_cuts = ensure_list(field_cuts)[:]

    def _get_list_of_grids(self):
        self._grids = self._base_region._grids

    def _is_fully_enclosed(self, grid):
        return False

    def _get_cut_mask(self, grid):
        point_mask = na.ones(grid.ActiveDimensions, dtype='bool')
        point_mask *= self._base_region._get_cut_mask(grid)
        for cut in self._field_cuts:
            point_mask *= eval(cut)
        return point_mask

class YTBooleanRegionBase(YTSelectionContainer3D):
    """
    A hybrid region built by boolean comparison between
    existing regions.
    """
    _type_name = "boolean"
    _con_args = ("regions")
    def __init__(self, regions, fields = None, pf = None, **kwargs):
        """
        This will build a hybrid region based on the boolean logic
        of the regions.

        Parameters
        ----------
        regions : list
            A list of region objects and strings describing the boolean logic
            to use when building the hybrid region. The boolean logic can be
            nested using parentheses.

        Examples
        --------
        >>> re1 = pf.h.region([0.5, 0.5, 0.5], [0.4, 0.4, 0.4],
            [0.6, 0.6, 0.6])
        >>> re2 = pf.h.region([0.5, 0.5, 0.5], [0.45, 0.45, 0.45],
            [0.55, 0.55, 0.55])
        >>> sp1 = pf.h.sphere([0.575, 0.575, 0.575], .03)
        >>> toroid_shape = pf.h.boolean([re1, "NOT", re2])
        >>> toroid_shape_with_hole = pf.h.boolean([re1, "NOT", "(", re2, "OR",
            sp1, ")"])
        """
        # Center is meaningless, but we'll define it all the same.
        YTSelectionContainer3D.__init__(self, [0.5]*3, fields, pf, **kwargs)
        self.regions = regions
        self._all_regions = []
        self._some_overlap = []
        self._all_overlap = []
        self._cut_masks = {}
        self._get_all_regions()
        self._make_overlaps()
        self._get_list_of_grids()

    def _get_all_regions(self):
        # Before anything, we simply find out which regions are involved in all
        # of this process, uniquely.
        for item in self.regions:
            if isinstance(item, types.StringType): continue
            self._all_regions.append(item)
            # So cut_masks don't get messed up.
            item._boolean_touched = True
        self._all_regions = na.unique(self._all_regions)

    def _make_overlaps(self):
        # Using the processed cut_masks, we'll figure out what grids
        # are left in the hybrid region.
        pbar = get_pbar("Building boolean", len(self._all_regions))
        for i, region in enumerate(self._all_regions):
            try:
                region._get_list_of_grids()
                alias = region
            except AttributeError:
                alias = region.data
            for grid in alias._grids:
                if grid in self._some_overlap or grid in self._all_overlap:
                    continue
                # Get the cut_mask for this grid in this region, and see
                # if there's any overlap with the overall cut_mask.
                overall = self._get_cut_mask(grid)
                local = force_array(alias._get_cut_mask(grid),
                    grid.ActiveDimensions)
                # Below we don't want to match empty masks.
                if overall.sum() == 0 and local.sum() == 0: continue
                # The whole grid is in the hybrid region if a) its cut_mask
                # in the original region is identical to the new one and b)
                # the original region cut_mask is all ones.
                if (local == na.bitwise_and(overall, local)).all() and \
                        (local == True).all():
                    self._all_overlap.append(grid)
                    continue
                if (overall == local).any():
                    # Some of local is in overall
                    self._some_overlap.append(grid)
                    continue
            pbar.update(i)
        pbar.finish()

    def __repr__(self):
        # We'll do this the slow way to be clear what's going on
        s = "%s (%s): " % (self.__class__.__name__, self.pf)
        s += "["
        for i, region in enumerate(self.regions):
            if region in ["OR", "AND", "NOT", "(", ")"]:
                s += region
            else:
                s += region.__repr__(clean = True)
            if i < (len(self.regions) - 1): s += ", "
        s += "]"
        return s

    def _is_fully_enclosed(self, grid):
        return (grid in self._all_overlap)

    def _get_list_of_grids(self):
        self._grids = na.array(self._some_overlap + self._all_overlap,
            dtype='object')

    def _get_cut_mask(self, grid, field=None):
        if self._is_fully_enclosed(grid):
            return True # We do not want child masking here
        if not isinstance(grid, (FakeGridForParticles,)) \
             and grid.id in self._cut_masks:
            return self._cut_masks[grid.id]
        # If we get this far, we have to generate the cut_mask.
        return self._get_level_mask(self.regions, grid)

    def _get_level_mask(self, ops, grid):
        level_masks = []
        end = 0
        for i, item in enumerate(ops):
            if end > 0 and i < end:
                # We skip over things inside parentheses on this level.
                continue
            if isinstance(item, YTDataContainer):
                # Add this regions cut_mask to level_masks
                level_masks.append(force_array(item._get_cut_mask(grid),
                    grid.ActiveDimensions))
            elif item == "AND" or item == "NOT" or item == "OR":
                level_masks.append(item)
            elif item == "(":
                # recurse down, and we'll append the results, which
                # should be a single cut_mask
                open_count = 0
                for ii, item in enumerate(ops[i + 1:]):
                    # We look for the matching closing parentheses to find
                    # where we slice ops.
                    if item == "(":
                        open_count += 1
                    if item == ")" and open_count > 0:
                        open_count -= 1
                    elif item == ")" and open_count == 0:
                        end = i + ii + 1
                        break
                level_masks.append(force_array(self._get_level_mask(ops[i + 1:end],
                    grid), grid.ActiveDimensions))
            elif isinstance(item.data, AMRData):
                level_masks.append(force_array(item.data._get_cut_mask(grid),
                    grid.ActiveDimensions))
            else:
                mylog.error("Item in the boolean construction unidentified.")
        # Now we do the logic on our level_mask.
        # There should be no nested logic anymore.
        # The first item should be a cut_mask,
        # so that will be our starting point.
        this_cut_mask = level_masks[0]
        for i, item in enumerate(level_masks):
            # I could use a slice above, but I'll keep i consistent instead.
            if i == 0: continue
            if item == "AND":
                # So, the next item in level_masks we want to AND.
                na.bitwise_and(this_cut_mask, level_masks[i+1], this_cut_mask)
            if item == "NOT":
                # It's convenient to remember that NOT == AND NOT
                na.bitwise_and(this_cut_mask, na.invert(level_masks[i+1]),
                    this_cut_mask)
            if item == "OR":
                na.bitwise_or(this_cut_mask, level_masks[i+1], this_cut_mask)
        if not isinstance(grid, FakeGridForParticles):
            self._cut_masks[grid.id] = this_cut_mask
        return this_cut_mask

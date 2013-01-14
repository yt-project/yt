"""
Data structures for Streaming, in-memory datasets

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

import weakref
import numpy as np

from yt.utilities.io_handler import io_registry
from yt.funcs import *
from yt.config import ytcfg
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridGeometryHandler
from yt.data_objects.static_output import \
    StaticOutput
from yt.utilities.logger import ytLogger as mylog
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.utilities.lib import \
    get_box_grids_level
from yt.utilities.decompose import \
    decompose_array, get_psize
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.flagging_methods import \
    FlaggingGrid

from .fields import \
    StreamFieldInfo, \
    add_stream_field, \
    KnownStreamFields

class StreamGrid(AMRGridPatch):
    """
    Class representing a single In-memory Grid instance.
    """

    __slots__ = ['proc_num']
    _id_offset = 0
    def __init__(self, id, hierarchy):
        """
        Returns an instance of StreamGrid with *id*, associated with *filename*
        and *hierarchy*.
        """
        #All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename = None, hierarchy = hierarchy)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    def _guess_properties_from_parent(self):
        rf = self.pf.refine_by
        my_ind = self.id - self._id_offset
        le = self.LeftEdge
        self.dds = self.Parent.dds/rf
        ParentLeftIndex = np.rint((self.LeftEdge-self.Parent.LeftEdge)/self.Parent.dds)
        self.start_index = rf*(ParentLeftIndex + self.Parent.get_global_startindex()).astype('int64')
        self.LeftEdge = self.Parent.LeftEdge + self.Parent.dds * ParentLeftIndex
        self.RightEdge = self.LeftEdge + self.ActiveDimensions*self.dds
        self.hierarchy.grid_left_edge[my_ind,:] = self.LeftEdge
        self.hierarchy.grid_right_edge[my_ind,:] = self.RightEdge
        self._child_mask = None
        self._child_index_mask = None
        self._child_indices = None
        self._setup_dx()

    def set_filename(self, filename):
        pass

    def __repr__(self):
        return "StreamGrid_%04i" % (self.id)

    @property
    def Parent(self):
        if self._parent_id == -1: return None
        return self.hierarchy.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.hierarchy.grids[cid - self._id_offset]
                for cid in self._children_ids]

class StreamHandler(object):
    def __init__(self, left_edges, right_edges, dimensions,
                 levels, parent_ids, particle_count, processor_ids,
                 fields, io = None, particle_types = {}):
        self.left_edges = left_edges
        self.right_edges = right_edges
        self.dimensions = dimensions
        self.levels = levels
        self.parent_ids = parent_ids
        self.particle_count = particle_count
        self.processor_ids = processor_ids
        self.num_grids = self.levels.size
        self.fields = fields
        self.io = io
        self.particle_types = particle_types
            
    def get_fields(self):
        return self.fields.all_fields

    def get_particle_type(self, field) :

        if self.particle_types.has_key(field) :
            return self.particle_types[field]
        else :
            return False
        
class StreamHierarchy(GridGeometryHandler):

    grid = StreamGrid

    def __init__(self, pf, data_style = None):
        self.data_style = data_style
        self.float_type = 'float64'
        self.parameter_file = weakref.proxy(pf) # for _obtain_enzo
        self.stream_handler = pf.stream_handler
        self.float_type = "float64"
        self.directory = os.getcwd()
        GridGeometryHandler.__init__(self, pf, data_style)

    def _count_grids(self):
        self.num_grids = self.stream_handler.num_grids

    def _parse_hierarchy(self):
        self.grid_dimensions = self.stream_handler.dimensions
        self.grid_left_edge[:] = self.stream_handler.left_edges
        self.grid_right_edge[:] = self.stream_handler.right_edges
        self.grid_levels[:] = self.stream_handler.levels
        self.grid_procs = self.stream_handler.processor_ids
        self.grid_particle_count[:] = self.stream_handler.particle_count
        mylog.debug("Copying reverse tree")
        self.grids = []
        # We enumerate, so it's 0-indexed id and 1-indexed pid
        for id in xrange(self.num_grids):
            self.grids.append(self.grid(id, self))
            self.grids[id].Level = self.grid_levels[id, 0]
        parent_ids = self.stream_handler.parent_ids
        if parent_ids is not None:
            reverse_tree = self.stream_handler.parent_ids.tolist()
            # Initial setup:
            for gid,pid in enumerate(reverse_tree):
                if pid >= 0:
                    self.grids[id]._parent_id = pid
                    self.grids[pid]._children_ids.append(self.grids[gid].id)
        else:
            mylog.debug("Reconstructing parent-child relationships")
            self._reconstruct_parent_child()
        self.max_level = self.grid_levels.max()
        mylog.debug("Preparing grids")
        temp_grids = np.empty(self.num_grids, dtype='object')
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = None
            grid._prepare_grid()
            grid.proc_num = self.grid_procs[i]
            temp_grids[i] = grid
        self.grids = temp_grids
        mylog.debug("Prepared")

    def _reconstruct_parent_child(self):
        mask = np.empty(len(self.grids), dtype='int32')
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[i,:],
                                self.grid_right_edge[i,:],
                                self.grid_levels[i] + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            ids = np.where(mask.astype("bool"))
            grid._children_ids = ids[0] # where is a tuple
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids): # Second pass
            for child in grid.Children:
                child._parent_id = i

    def _initialize_grid_arrays(self):
        GridGeometryHandler._initialize_grid_arrays(self)
        self.grid_procs = np.zeros((self.num_grids,1),'int32')

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridGeometryHandler._setup_classes(self, dd)

    def _detect_fields(self):
        self.field_list = list(set(self.stream_handler.get_fields()))

    def _populate_grid_objects(self):
        for g in self.grids:
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _setup_data_io(self):
        if self.stream_handler.io is not None:
            self.io = self.stream_handler.io
        else:
            self.io = io_registry[self.data_style](self.stream_handler)

    def update_data(self, data) :

        """
        Update the stream data with a new data dict. If fields already exist,
        they will be replaced, but if they do not, they will be added. Fields
        already in the stream but not part of the data dict will be left
        alone. 
        """
        
        particle_types = set_particle_types(data[0])

        for key in data[0].keys() :
            if key is "number_of_particles": continue
            self.stream_handler.particle_types[key] = particle_types[key]
            if key not in self.field_list:
                self.field_list.append(key)
                
        for i, grid in enumerate(self.grids) :
            if data[i].has_key("number_of_particles") :
                grid.NumberOfParticles = data[i].pop("number_of_particles")
            for key in data[i].keys() :
                if key in grid.keys() : grid.field_data.pop(key, None)
                self.stream_handler.fields[grid.id][key] = data[i][key]
            
        self._detect_fields()
        self._setup_unknown_fields()
                
class StreamStaticOutput(StaticOutput):
    _hierarchy_class = StreamHierarchy
    _fieldinfo_fallback = StreamFieldInfo
    _fieldinfo_known = KnownStreamFields
    _data_style = 'stream'

    def __init__(self, stream_handler, storage_filename = None):
        #if parameter_override is None: parameter_override = {}
        #self._parameter_override = parameter_override
        #if conversion_override is None: conversion_override = {}
        #self._conversion_override = conversion_override

        self.stream_handler = stream_handler
        StaticOutput.__init__(self, "InMemoryParameterFile", self._data_style)
        self.storage_filename = storage_filename

        self.units = {}
        self.time_units = {}

    def _parse_parameter_file(self):
        self.basename = self.stream_handler.name
        self.parameters['CurrentTimeIdentifier'] = time.time()
        self.unique_identifier = self.parameters["CurrentTimeIdentifier"]
        self.domain_left_edge = self.stream_handler.domain_left_edge[:]
        self.domain_right_edge = self.stream_handler.domain_right_edge[:]
        self.refine_by = self.stream_handler.refine_by
        self.dimensionality = self.stream_handler.dimensionality
        self.domain_dimensions = self.stream_handler.domain_dimensions
        self.current_time = self.stream_handler.simulation_time
        if self.stream_handler.cosmology_simulation:
            self.cosmological_simulation = 1
            self.current_redshift = self.stream_handler.current_redshift
            self.omega_lambda = self.stream_handler.omega_lambda
            self.omega_matter = self.stream_handler.omega_matter
            self.hubble_constant = self.stream_handler.hubble_constant
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _set_units(self):
        pass

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return False

class StreamDictFieldHandler(dict):

    @property
    def all_fields(self): return self[0].keys()

def set_particle_types(data) :

    particle_types = {}
    
    for key in data.keys() :

        if key is "number_of_particles": continue
        
        if len(data[key].shape) == 1:
            particle_types[key] = True
        else :
            particle_types[key] = False
    
    return particle_types

def assign_particle_data(pf, pdata) :

    """
    Assign particle data to the grids using find_points. This
    will overwrite any existing particle data, so be careful!
    """
    
    if pf.h.num_grids > 1 :

        try:
            x, y, z = (pdata["all","particle_position_%s" % ax] for ax in 'xyz')
        except KeyError:
            raise KeyError("Cannot decompose particle data without position fields!")
        
        particle_grids, particle_grid_inds = pf.h.find_points(x,y,z)
        idxs = np.argsort(particle_grid_inds)
        particle_grid_count = np.bincount(particle_grid_inds,
                                          minlength=pf.h.num_grids)
        particle_indices = np.zeros(pf.h.num_grids + 1, dtype='int64')
        if pf.h.num_grids > 1 :
            np.add.accumulate(particle_grid_count.squeeze(),
                              out=particle_indices[1:])
        else :
            particle_indices[1] = particle_grid_count.squeeze()
    
        pdata.pop("number_of_particles")    
        grid_pdata = []
        
        for i, pcount in enumerate(particle_grid_count) :
            grid = {}
            grid["number_of_particles"] = pcount
            start = particle_indices[i]
            end = particle_indices[i+1]
            for key in pdata.keys() :
                grid[key] = pdata[key][idxs][start:end]
            grid_pdata.append(grid)

    else :

        grid_pdata = [pdata]
        
    pf.h.update_data(grid_pdata)
                                        
def load_uniform_grid(data, domain_dimensions, sim_unit_to_cm, bbox=None,
                      nprocs=1, sim_time=0.0):
    r"""Load a uniform grid of data into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamHandler`.

    This should allow a uniform grid of data to be loaded directly into yt and
    analyzed as would any others.  This comes with several caveats:
        * Units will be incorrect unless the data has already been converted to
          cgs.
        * Some functions may behave oddly, and parallelism will be
          disappointing or non-existent in most cases.
        * Particles may be difficult to integrate.

    Particle fields are detected as one-dimensional fields. The number of particles
    is set by the "number_of_particles" key in data.
    
    Parameters
    ----------
    data : dict
        This is a dict of numpy arrays, where the keys are the field names.
    domain_dimensions : array_like
        This is the domain dimensions of the grid
    sim_unit_to_cm : float
        Conversion factor from simulation units to centimeters
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units sim_unit_to_cm
    nprocs: integer, optional
        If greater than 1, will create this number of subarrays out of data
    sim_time : float, optional
        The simulation time in seconds

    Examples
    --------

    >>> arr = np.random.random((128, 128, 129))
    >>> data = dict(Density = arr)
    >>> bbox = np.array([[0., 1.0], [-1.5, 1.5], [1.0, 2.5]])
    >>> pf = load_uniform_grid(data, arr.shape, 3.08e24, bbox=bbox, nprocs=12)

    """

    domain_dimensions = np.array(domain_dimensions)
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], 'float64')
    domain_left_edge = np.array(bbox[:, 0], 'float64')
    domain_right_edge = np.array(bbox[:, 1], 'float64')
    grid_levels = np.zeros(nprocs, dtype='int32').reshape((nprocs,1))

    sfh = StreamDictFieldHandler()
    
    if data.has_key("number_of_particles") :
        number_of_particles = data.pop("number_of_particles")
    else :
        number_of_particles = int(0)
    
    if number_of_particles > 0 :
        particle_types = set_particle_types(data)
        pdata = {}
        pdata["number_of_particles"] = number_of_particles
        for key in data.keys() :
            if len(data[key].shape) == 1 :
                pdata[key] = data.pop(key)
    else :
        particle_types = {}
    
    if nprocs > 1:
        temp = {}
        new_data = {}
        for key in data.keys():
            psize = get_psize(np.array(data[key].shape), nprocs)
            grid_left_edges, grid_right_edges, temp[key] = \
                             decompose_array(data[key], psize, bbox)
            grid_dimensions = np.array([grid.shape for grid in temp[key]],
                                       dtype="int32")
        for gid in range(nprocs):
            new_data[gid] = {}
            for key in temp.keys():
                new_data[gid].update({key:temp[key][gid]})
        sfh.update(new_data)
        del new_data, temp
    else:
        sfh.update({0:data})
        grid_left_edges = domain_left_edge
        grid_right_edges = domain_right_edge
        grid_dimensions = domain_dimensions.reshape(nprocs,3).astype("int32")

    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        -np.ones(nprocs, dtype='int64'),
        np.zeros(nprocs, dtype='int64').reshape(nprocs,1), # Temporary
        np.zeros(nprocs).reshape((nprocs,1)),
        sfh,
        particle_types=particle_types
    )

    handler.name = "UniformGridData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    handler.dimensionality = 3
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    spf = StreamStaticOutput(handler)
    spf.units["cm"] = sim_unit_to_cm
    spf.units['1'] = 1.0
    spf.units["unitary"] = 1.0
    box_in_mpc = sim_unit_to_cm / mpc_conversion['cm']
    for unit in mpc_conversion.keys():
        spf.units[unit] = mpc_conversion[unit] * box_in_mpc

    # Now figure out where the particles go

    if number_of_particles > 0 :
        if ("all", "particle_position_x") not in pdata:
            pdata_ftype = {}
            for f in [k for k in sorted(pdata)]:
                if not hasattr(pdata[f], "shape"): continue
                mylog.debug("Reassigning '%s' to ('all','%s')", f, f)
                pdata_ftype["all",f] = pdata.pop(f)
            pdata_ftype.update(pdata)
            pdata = pdata_ftype
        assign_particle_data(spf, pdata)
    
    return spf

def load_amr_grids(grid_data, domain_dimensions, sim_unit_to_cm, bbox=None,
                   sim_time=0.0):
    r"""Load a set of grids of data into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamHandler`.

    This should allow a sequence of grids of varying resolution of data to be
    loaded directly into yt and analyzed as would any others.  This comes with
    several caveats:
        * Units will be incorrect unless the data has already been converted to
          cgs.
        * Some functions may behave oddly, and parallelism will be
          disappointing or non-existent in most cases.
        * Particles may be difficult to integrate.
        * No consistency checks are performed on the hierarchy

    Parameters
    ----------
    grid_data : list of dicts
        This is a list of dicts.  Each dict must have entries "left_edge",
        "right_edge", "dimensions", "level", and then any remaining entries are
        assumed to be fields.  They also may include a particle count, otherwise
        assumed to be zero. This will be modified in place and can't be
        assumed to be static.
    domain_dimensions : array_like
        This is the domain dimensions of the grid
    sim_unit_to_cm : float
        Conversion factor from simulation units to centimeters
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units sim_unit_to_cm
    sim_time : float, optional
        The simulation time in seconds

    Examples
    --------

    >>> grid_data = [
    ...     dict(left_edge = [0.0, 0.0, 0.0],
    ...          right_edge = [1.0, 1.0, 1.],
    ...          level = 0,
    ...          dimensions = [32, 32, 32],
    ...          number_of_particles = 0)
    ...     dict(left_edge = [0.25, 0.25, 0.25],
    ...          right_edge = [0.75, 0.75, 0.75],
    ...          level = 1,
    ...          dimensions = [32, 32, 32],
    ...          number_of_particles = 0)
    ... ]
    ... 
    >>> for g in grid_data:
    ...     g["Density"] = np.random.random(g["dimensions"]) * 2**g["level"]
    ...
    >>> pf = load_amr_grids(grid_data, [32, 32, 32], 1.0)
    """

    domain_dimensions = np.array(domain_dimensions)
    ngrids = len(grid_data)
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], 'float64')
    domain_left_edge = np.array(bbox[:, 0], 'float64')
    domain_right_edge = np.array(bbox[:, 1], 'float64')
    grid_levels = np.zeros((ngrids, 1), dtype='int32')
    grid_left_edges = np.zeros((ngrids, 3), dtype="float32")
    grid_right_edges = np.zeros((ngrids, 3), dtype="float32")
    grid_dimensions = np.zeros((ngrids, 3), dtype="int32")
    number_of_particles = np.zeros((ngrids,1), dtype='int64')
    sfh = StreamDictFieldHandler()
    for i, g in enumerate(grid_data):
        grid_left_edges[i,:] = g.pop("left_edge")
        grid_right_edges[i,:] = g.pop("right_edge")
        grid_dimensions[i,:] = g.pop("dimensions")
        grid_levels[i,:] = g.pop("level")
        if g.has_key("number_of_particles") :
            number_of_particles[i,:] = g.pop("number_of_particles")  
        sfh[i] = g
            
    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        None, # parent_ids is none
        number_of_particles,
        np.zeros(ngrids).reshape((ngrids,1)),
        sfh,
        particle_types=set_particle_types(grid_data[0])
    )

    handler.name = "AMRGridData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    handler.dimensionality = 3
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    spf = StreamStaticOutput(handler)
    spf.units["cm"] = sim_unit_to_cm
    spf.units['1'] = 1.0
    spf.units["unitary"] = 1.0
    box_in_mpc = sim_unit_to_cm / mpc_conversion['cm']
    for unit in mpc_conversion.keys():
        spf.units[unit] = mpc_conversion[unit] * box_in_mpc
    return spf

def refine_amr(base_pf, refinement_criteria, fluid_operators, max_level,
               callback = None):
    r"""Given a base parameter file, repeatedly apply refinement criteria and
    fluid operators until a maximum level is reached.

    Parameters
    ----------
    base_pf : StaticOutput
        This is any static output.  It can also be a stream static output, for
        instance as returned by load_uniform_data.
    refinement_critera : list of :class:`~yt.utilities.flagging_methods.FlaggingMethod`
        These criteria will be applied in sequence to identify cells that need
        to be refined.
    fluid_operators : list of :class:`~yt.utilities.initial_conditions.FluidOperator`
        These fluid operators will be applied in sequence to all resulting
        grids.
    max_level : int
        The maximum level to which the data will be refined
    callback : function, optional
        A function that will be called at the beginning of each refinement
        cycle, with the current parameter file.

    Examples
    --------
    >>> domain_dims = (32, 32, 32)
    >>> data = np.zeros(domain_dims) + 0.25
    >>> fo = [ic.CoredSphere(0.05, 0.3, [0.7,0.4,0.75], {"Density": (0.25, 100.0)})]
    >>> rc = [fm.flagging_method_registry["overdensity"](8.0)]
    >>> ug = load_uniform_grid({'Density': data}, domain_dims, 1.0)
    >>> pf = refine_amr(ug, rc, fo, 5)
    """

    # If we have particle data, set it aside for now

    number_of_particles = np.sum([grid.NumberOfParticles
                                  for grid in base_pf.h.grids])

    if number_of_particles > 0 :
        pdata = {}
        for field in base_pf.h.field_list :
            if base_pf.field_info[field].particle_type :
                pdata[field] = np.concatenate([grid[field]
                                               for grid in base_pf.h.grids])
        pdata["number_of_particles"] = number_of_particles
        
    last_gc = base_pf.h.num_grids
    cur_gc = -1
    pf = base_pf    
    bbox = np.array( [ (pf.domain_left_edge[i], pf.domain_right_edge[i])
                       for i in range(3) ])
    while pf.h.max_level < max_level and last_gc != cur_gc:
        mylog.info("Refining another level.  Current max level: %s",
                  pf.h.max_level)
        last_gc = pf.h.grids.size
        for m in fluid_operators: m.apply(pf)
        if callback is not None: callback(pf)
        grid_data = []
        for g in pf.h.grids:
            gd = dict( left_edge = g.LeftEdge,
                       right_edge = g.RightEdge,
                       level = g.Level,
                       dimensions = g.ActiveDimensions )
            for field in pf.h.field_list:
                if not pf.field_info[field].particle_type :
                    gd[field] = g[field]
            grid_data.append(gd)
            if g.Level < pf.h.max_level: continue
            fg = FlaggingGrid(g, refinement_criteria)
            nsg = fg.find_subgrids()
            for sg in nsg:
                LE = sg.left_index * g.dds + pf.domain_left_edge
                dims = sg.dimensions * pf.refine_by
                grid = pf.h.smoothed_covering_grid(g.Level + 1, LE, dims)
                gd = dict(left_edge = LE, right_edge = grid.right_edge,
                          level = g.Level + 1, dimensions = dims)
                for field in pf.h.field_list:
                    if not pf.field_info[field].particle_type :
                        gd[field] = grid[field]
                grid_data.append(gd)
        
        pf = load_amr_grids(grid_data, pf.domain_dimensions, 1.0,
                            bbox = bbox)
        cur_gc = pf.h.num_grids

    # Now reassign particle data to grids

    if number_of_particles > 0:
        if ("all", "particle_position_x") not in pdata:
            pdata_ftype = {}
            for f in [k for k in sorted(pdata)]:
                if not hasattr(pdata[f], "shape"): continue
                mylog.debug("Reassigning '%s' to ('all','%s')", f, f)
                pdata_ftype["all",f] = pdata.pop(f)
            pdata_ftype.update(pdata)
            pdata = pdata_ftype
        assign_particle_data(pf, pdata)
    
    return pf

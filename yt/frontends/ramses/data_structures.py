"""
RAMSES-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as na
import stat
import weakref

from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.data_objects.hierarchy import \
      AMRHierarchy
from yt.data_objects.static_output import \
      StaticOutput

try:
    import _ramses_reader
except ImportError:
    _ramses_reader = None
from .fields import RAMSESFieldInfo, KnownRAMSESFields
from yt.utilities.definitions import \
    mpc_conversion
from yt.utilities.amr_utils import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

def num_deep_inc(f):
    def wrap(self, *args, **kwargs):
        self.num_deep += 1
        rv = f(self, *args, **kwargs)
        self.num_deep -= 1
        return rv
    return wrap

class RAMSESGrid(AMRGridPatch):
    _id_offset = 0
    #__slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, hierarchy, level, locations, start_index):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Level = level
        self.Parent = []
        self.Children = []
        self.locations = locations
        self.domain = locations[0,0]
        self.start_index = start_index.copy()

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.
        """
        if self.start_index != None:
            return self.start_index
        if len(self.Parent) == 0:
            start_index = self.LeftEdge / self.dds
            return na.rint(start_index).astype('int64').ravel()
        pdx = self.Parent[0].dds
        start_index = (self.Parent[0].get_global_startindex()) + \
                       na.rint((self.LeftEdge - self.Parent[0].LeftEdge)/pdx)
        self.start_index = (start_index*self.pf.refine_by).astype('int64').ravel()
        return self.start_index

    def __repr__(self):
        return "RAMSESGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class RAMSESHierarchy(AMRHierarchy):

    grid = RAMSESGrid
    _handle = None
    
    def __init__(self,pf,data_style='ramses'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.tree_proxy = pf.ramses_tree

        self.float_type = na.float64
        AMRHierarchy.__init__(self,pf,data_style)

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        self.field_list = self.tree_proxy.field_names[:]
    
    def _setup_field_list(self):
        if self.parameter_file.use_particles:
            # We know which particle fields will exist -- pending further
            # changes in the future.
            for field in art_particle_field_names:
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
                # Note that we call add_field on the field_info directly.  This
                # will allow the same field detection mechanism to work for 1D,
                # 2D and 3D fields.
                self.pf.field_info.add_field(field, NullFunc,
                                             convert_function=cf,
                                             take_log=False, particle_type=True)

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        # We have to do all the patch-coalescing here.
        LEVEL_OF_EDGE = 7
        MAX_EDGE = (2 << (LEVEL_OF_EDGE- 1))
        level_info = self.tree_proxy.count_zones()
        num_ogrids = sum(level_info)
        ogrid_left_edge = na.zeros((num_ogrids,3), dtype='float64')
        ogrid_right_edge = na.zeros((num_ogrids,3), dtype='float64')
        ogrid_levels = na.zeros((num_ogrids,1), dtype='int32')
        ogrid_file_locations = na.zeros((num_ogrids,6), dtype='int64')
        ogrid_hilbert_indices = na.zeros(num_ogrids, dtype='uint64')
        ochild_masks = na.zeros((num_ogrids, 8), dtype='int32')
        self.tree_proxy.fill_hierarchy_arrays(
            self.pf.domain_dimensions,
            ogrid_left_edge, ogrid_right_edge,
            ogrid_levels, ogrid_file_locations, 
            ogrid_hilbert_indices, ochild_masks,
            self.pf.domain_left_edge, self.pf.domain_right_edge)
        #insert_ipython()
        # Now we can rescale
        mi, ma = ogrid_left_edge.min(), ogrid_right_edge.max()
        DL = self.pf.domain_left_edge
        DR = self.pf.domain_right_edge
        DW = DR - DL
        ogrid_left_edge = (ogrid_left_edge - mi)/(ma - mi) * (DR - DL) + DL
        ogrid_right_edge = (ogrid_right_edge - mi)/(ma - mi) * (DR - DL) + DL
        #import pdb;pdb.set_trace()
        # We now have enough information to run the patch coalescing 
        self.proto_grids = []
        for level in xrange(len(level_info)):
            if level_info[level] == 0: continue
            # Get the indices of grids on this level
            ggi = (ogrid_levels == level).ravel()
            dims = na.ones((ggi.sum(), 3), dtype='int64') * 2 
            mylog.info("Re-gridding level %s: %s octree grids", level, ggi.sum())
            nd = self.pf.domain_dimensions * 2**level
            fl = ogrid_file_locations[ggi,:]
            # Now our initial protosubgrid
            #if level == 6: raise RuntimeError
            # We want grids that cover no more than MAX_EDGE cells in every direction
            psgs = []
            # left_index is integers of the index, with respect to this level
            left_index = na.rint((ogrid_left_edge[ggi,:]) * nd / DW ).astype('int64')
            # we've got octs, so it's +2
            pbar = get_pbar("Re-gridding ", left_index.shape[0])
            dlp = [None, None, None]
            i = 0
            # We now calculate the hilbert curve position of every left_index,
            # of the octs, with respect to a lower order hilbert curve.
            left_index_gridpatch = left_index >> LEVEL_OF_EDGE
            order = max(level + 1 - LEVEL_OF_EDGE, 0)
            # I'm not sure the best way to do this.
            hilbert_indices = _ramses_reader.get_hilbert_indices(order, left_index_gridpatch)
            #print level, hilbert_indices.min(), hilbert_indices.max()
            # Strictly speaking, we don't care about the index of any
            # individual oct at this point.  So we can then split them up.
            unique_indices = na.unique(hilbert_indices)
            mylog.debug("Level % 2i has % 10i unique indices for %0.3e octs",
                        level, unique_indices.size, hilbert_indices.size)
            locs, lefts = _ramses_reader.get_array_indices_lists(
                        hilbert_indices, unique_indices, left_index, fl)
            for ddleft_index, ddfl in zip(lefts, locs):
                for idomain in na.unique(ddfl[:,0]):
                    dom_ind = ddfl[:,0] == idomain
                    dleft_index = ddleft_index[dom_ind,:]
                    dfl = ddfl[dom_ind,:]
                    initial_left = na.min(dleft_index, axis=0)
                    idims = (na.max(dleft_index, axis=0) - initial_left).ravel()+2
                    psg = _ramses_reader.ProtoSubgrid(initial_left, idims,
                                    dleft_index, dfl)
                    if psg.efficiency <= 0: continue
                    self.num_deep = 0
                    psgs.extend(_ramses_reader.recursive_patch_splitting(
                        psg, idims, initial_left, 
                        dleft_index, dfl))
            mylog.debug("Done with level % 2i", level)
            pbar.finish()
            self.proto_grids.append(psgs)
            print sum(len(psg.grid_file_locations) for psg in psgs)
            sums = na.zeros(3, dtype='int64')
            mylog.info("Final grid count: %s", len(self.proto_grids[level]))
            if len(self.proto_grids[level]) == 1: continue
            #for g in self.proto_grids[level]:
            #    sums += [s.sum() for s in g.sigs]
            #assert(na.all(sums == dims.prod(axis=1).sum()))
        self.num_grids = sum(len(l) for l in self.proto_grids)

    def _parse_hierarchy(self):
        # We have important work to do
        grids = []
        gi = 0
        DL, DR = self.pf.domain_left_edge, self.pf.domain_right_edge
        DW = DR - DL
        for level, grid_list in enumerate(self.proto_grids):
            for g in grid_list:
                fl = g.grid_file_locations
                props = g.get_properties()
                self.grid_left_edge[gi,:] = (props[0,:] / (2.0**(level+1))) * DW + DL
                self.grid_right_edge[gi,:] = (props[1,:] / (2.0**(level+1))) * DW + DL
                self.grid_dimensions[gi,:] = props[2,:]
                self.grid_levels[gi,:] = level
                grids.append(self.grid(gi, self, level, fl, props[0,:]))
                gi += 1
        self.proto_grids = []
        self.grids = na.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self.grids[gi] = g

    def _populate_grid_objects(self):
        mask = na.empty(self.grids.size, dtype='int32')
        print self.grid_levels.dtype
        for gi,g in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[gi,:],
                                self.grid_right_edge[gi,:],
                                g.Level - 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            parents = self.grids[mask.astype("bool")]
            if len(parents) > 0:
                g.Parent.extend((p for p in parents.tolist()
                        if p.locations[0,0] == g.locations[0,0]))
                for p in parents: p.Children.append(g)
            # Now we do overlapping siblings; note that one has to "win" with
            # siblings, so we assume the lower ID one will "win"
            get_box_grids_level(self.grid_left_edge[gi,:],
                                self.grid_right_edge[gi,:],
                                g.Level,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask, gi)
            mask[gi] = False
            siblings = self.grids[mask.astype("bool")]
            if len(siblings) > 0:
                g.OverlappingSiblings = siblings.tolist()
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.tree_proxy)

class RAMSESStaticOutput(StaticOutput):
    _hierarchy_class = RAMSESHierarchy
    _fieldinfo_fallback = RAMSESFieldInfo
    _fieldinfo_known = KnownRAMSESFields
    _handle = None
    
    def __init__(self, filename, data_style='ramses',
                 storage_filename = None):
        # Here we want to initiate a traceback, if the reader is not built.
        import _ramses_reader
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
        
    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = self.parameters['unit_t']
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        self.time_units['Myr'] = self.time_units['years'] / 1.0e6
        self.time_units['Gyr']  = self.time_units['years'] / 1.0e9
        self.conversion_factors["Density"] = self.parameters['unit_d']
        vel_u = self.parameters['unit_l'] / self.parameters['unit_t']
        self.conversion_factors["x-velocity"] = vel_u
        self.conversion_factors["y-velocity"] = vel_u
        self.conversion_factors["z-velocity"] = vel_u

    def _setup_nounits_units(self):
        for unit in mpc_conversion.keys():
            self.units[unit] = self.parameters['unit_l'] * mpc_conversion[unit] / mpc_conversion["cm"]

    def _parse_parameter_file(self):
        # hardcoded for now
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'ramses'
        self.parameters["Time"] = 1. # default unit is 1...

        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self.ramses_tree = _ramses_reader.RAMSES_tree_proxy(self.parameter_filename)
        rheader = self.ramses_tree.get_file_info()
        self.parameters.update(rheader)
        self.current_time = self.parameters['time'] * self.parameters['unit_t']
        self.domain_right_edge = na.ones(3, dtype='float64') \
                                           * rheader['boxlen']
        self.domain_left_edge = na.zeros(3, dtype='float64')
        self.domain_dimensions = na.ones(3, dtype='int32') * 2
        # This is likely not true, but I am not sure how to otherwise
        # distinguish them.
        mylog.warning("No current mechanism of distinguishing cosmological simulations in RAMSES!")
        self.cosmological_simulation = 1
        self.current_redshift = (1.0 / rheader["aexp"]) - 1.0
        self.omega_lambda = rheader["omega_l"]
        self.omega_matter = rheader["omega_m"]
        self.hubble_constant = rheader["H0"]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not os.path.basename(args[0]).startswith("info_"): return False
        fn = args[0].replace("info_", "amr_").replace(".txt", ".out00001")
        print fn
        return os.path.exists(fn)


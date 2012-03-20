"""
Data structures for Enzo

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
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

import h5py
import weakref
import numpy as na
import os
import stat
import string
import re
try:
    from pyhdf_np import SD
except ImportError:
    pass

from itertools import izip

from yt.funcs import *
from yt.config import ytcfg
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.hierarchy import \
    AMRHierarchy
from yt.data_objects.static_output import \
    StaticOutput
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.utilities.definitions import mpc_conversion
from yt.utilities import hdf5_light_reader
from yt.utilities.logger import ytLogger as mylog

from .definitions import parameterDict
from .fields import \
    EnzoFieldInfo, Enzo2DFieldInfo, Enzo1DFieldInfo, \
    add_enzo_field, add_enzo_2d_field, add_enzo_1d_field, \
    KnownEnzoFields

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_blocking_call

class EnzoGrid(AMRGridPatch):
    """
    Class representing a single Enzo Grid instance.
    """

    __slots__ = []
    def __init__(self, id, hierarchy):
        """
        Returns an instance of EnzoGrid with *id*, associated with
        *filename* and *hierarchy*.
        """
        #All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename = None, hierarchy = hierarchy)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    def _guess_properties_from_parent(self):
        """
        We know that our grid boundary occurs on the cell boundary of our
        parent.  This can be a very expensive process, but it is necessary
        in some hierarchys, where yt is unable to generate a completely
        space-filling tiling of grids, possibly due to the finite accuracy in a
        standard Enzo hierarchy file.
        """
        rf = self.pf.refine_by
        my_ind = self.id - self._id_offset
        le = self.LeftEdge
        self.dds = self.Parent.dds/rf
        ParentLeftIndex = na.rint((self.LeftEdge-self.Parent.LeftEdge)/self.Parent.dds)
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
        """
        Intelligently set the filename.
        """
        if self.hierarchy._strip_path:
            self.filename = os.path.join(self.hierarchy.directory,
                                         os.path.basename(filename))
        elif filename[0] == os.path.sep:
            self.filename = filename
        else:
            self.filename = os.path.join(self.hierarchy.directory, filename)
        return

    def __repr__(self):
        return "EnzoGrid_%04i" % (self.id)

    @property
    def Parent(self):
        if self._parent_id == -1: return None
        return self.hierarchy.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.hierarchy.grids[cid - self._id_offset]
                for cid in self._children_ids]

class EnzoGridInMemory(EnzoGrid):
    __slots__ = ['proc_num']
    def set_filename(self, filename):
        pass

class EnzoGridGZ(EnzoGrid):

    __slots__ = ()

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False,
                             smoothed=False):
        # We ignore smoothed in this case.
        if n_zones > 3:
            return EnzoGrid.retrieve_ghost_zones(
                self, n_zones, fields, all_levels, smoothed)
        # ----- Below is mostly the original code, except we remove the field
        # ----- access section
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        nl = self.get_global_startindex() - n_zones
        nr = nl + self.ActiveDimensions + 2*n_zones
        new_left_edge = nl * self.dds + self.pf.domain_left_edge
        new_right_edge = nr * self.dds + self.pf.domain_left_edge
        # Something different needs to be done for the root grid, though
        level = self.Level
        args = (level, new_left_edge, new_right_edge)
        kwargs = {'dims': self.ActiveDimensions + 2*n_zones,
                  'num_ghost_zones':n_zones,
                  'use_pbar':False}
        # This should update the arguments to set the field parameters to be
        # those of this grid.
        kwargs.update(self.field_parameters)
        if smoothed:
            #cube = self.hierarchy.smoothed_covering_grid(
            #    level, new_left_edge, new_right_edge, **kwargs)
            cube = self.hierarchy.smoothed_covering_grid(
                level, new_left_edge, **kwargs)
        else:
            cube = self.hierarchy.covering_grid(
                level, new_left_edge, **kwargs)
        # ----- This is EnzoGrid.get_data, duplicated here mostly for
        # ----  efficiency's sake.
        sl = [slice(3 - n_zones, -(3 - n_zones)) for i in range(3)]
        if fields is None: return cube
        for field in ensure_list(fields):
            if field in self.hierarchy.field_list:
                conv_factor = 1.0
                if self.pf.field_info.has_key(field):
                    conv_factor = self.pf.field_info[field]._convert_function(self)
                if self.pf.field_info[field].particle_type: continue
                temp = self.hierarchy.io._read_raw_data_set(self, field)
                temp = temp.swapaxes(0, 2)
                cube.field_data[field] = na.multiply(temp, conv_factor, temp)[sl]
        return cube

class EnzoHierarchy(AMRHierarchy):

    _strip_path = False
    grid = EnzoGrid

    def __init__(self, pf, data_style):
        
        self.data_style = data_style
        if pf.file_style != None:
            self._bn = pf.file_style
        else:
            self._bn = "%s.cpu%%04i"
        self.hierarchy_filename = os.path.abspath(
            "%s.hierarchy" % (pf.parameter_filename))
        harray_fn = self.hierarchy_filename[:-9] + "harrays"
        if ytcfg.getboolean("yt","serialize") and os.path.exists(harray_fn):
            try:
                harray_fp = h5py.File(harray_fn)
                self.num_grids = harray_fp["/Level"].len()
                harray_fp.close()
            except IOError:
                pass
        elif os.path.getsize(self.hierarchy_filename) == 0:
            raise IOError(-1,"File empty", self.hierarchy_filename)
        self.directory = os.path.dirname(self.hierarchy_filename)

        # For some reason, r8 seems to want Float64
        if pf.has_key("CompilerPrecision") \
            and pf["CompilerPrecision"] == "r4":
            self.float_type = 'float32'
        else:
            self.float_type = 'float64'

        AMRHierarchy.__init__(self, pf, data_style)
        # sync it back
        self.parameter_file.data_style = self.data_style

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        test_grid = test_grid_id = None
        self.num_stars = 0
        for line in rlines(open(self.hierarchy_filename, "rb")):
            if line.startswith("BaryonFileName") or \
               line.startswith("ParticleFileName") or \
               line.startswith("FileName "):
                test_grid = line.split("=")[-1].strip().rstrip()
            if line.startswith("NumberOfStarParticles"):
                self.num_stars = int(line.split("=")[-1])
            if line.startswith("Grid "):
                self.num_grids = test_grid_id = int(line.split("=")[-1])
                break
        self._guess_data_style(self.pf.dimensionality, test_grid, test_grid_id)

    def _guess_data_style(self, rank, test_grid, test_grid_id):
        if test_grid[0] != os.path.sep:
            test_grid = os.path.join(self.directory, test_grid)
        if not os.path.exists(test_grid):
            test_grid = os.path.join(self.directory,
                                    os.path.basename(test_grid))
            mylog.debug("Your data uses the annoying hardcoded path.")
            self._strip_path = True
        if self.data_style is not None: return
        try:
            a = SD.SD(test_grid)
            self.data_style = 'enzo_hdf4'
            mylog.debug("Detected HDF4")
        except:
            try:
                list_of_sets = hdf5_light_reader.ReadListOfDatasets(test_grid, "/")
            except:
                print "Could not find dataset.  Defaulting to packed HDF5"
                list_of_sets = []
            if len(list_of_sets) == 0 and rank == 3:
                mylog.debug("Detected packed HDF5")
                if self.parameters.get("WriteGhostZones", 0) == 1:
                    self.data_style= "enzo_packed_3d_gz"
                    self.grid = EnzoGridGZ
                else:
                    self.data_style = 'enzo_packed_3d'
            elif len(list_of_sets) > 0 and rank == 3:
                mylog.debug("Detected unpacked HDF5")
                self.data_style = 'enzo_hdf5'
            elif len(list_of_sets) == 0 and rank == 2:
                mylog.debug("Detect packed 2D")
                self.data_style = 'enzo_packed_2d'
            elif len(list_of_sets) == 0 and rank == 1:
                mylog.debug("Detect packed 1D")
                self.data_style = 'enzo_packed_1d'
            else:
                raise TypeError

    # Sets are sorted, so that won't work!
    def _parse_hierarchy(self):
        def _next_token_line(token, f):
            for line in f:
                if line.startswith(token):
                    return line.split()[2:]
        if os.path.exists(self.hierarchy_filename[:-9] + "harrays"):
            if self._parse_binary_hierarchy(): return
        t1 = time.time()
        pattern = r"Pointer: Grid\[(\d*)\]->NextGrid(Next|This)Level = (\d*)\s+$"
        patt = re.compile(pattern)
        f = open(self.hierarchy_filename, "rb")
        self.grids = [self.grid(1, self)]
        self.grids[0].Level = 0
        si, ei, LE, RE, fn, np = [], [], [], [], [], []
        all = [si, ei, LE, RE, fn]
        pbar = get_pbar("Parsing Hierarchy", self.num_grids)
        for grid_id in xrange(self.num_grids):
            pbar.update(grid_id)
            # We will unroll this list
            si.append(_next_token_line("GridStartIndex", f))
            ei.append(_next_token_line("GridEndIndex", f))
            LE.append(_next_token_line("GridLeftEdge", f))
            RE.append(_next_token_line("GridRightEdge", f))
            nb = int(_next_token_line("NumberOfBaryonFields", f)[0])
            fn.append(["-1"])
            if nb > 0: fn[-1] = _next_token_line("BaryonFileName", f)
            np.append(int(_next_token_line("NumberOfParticles", f)[0]))
            if nb == 0 and np[-1] > 0: fn[-1] = _next_token_line("ParticleFileName", f)
            for line in f:
                if len(line) < 2: break
                if line.startswith("Pointer:"):
                    vv = patt.findall(line)[0]
                    self.__pointer_handler(vv)
        pbar.finish()
        self._fill_arrays(ei, si, LE, RE, np)
        temp_grids = na.empty(self.num_grids, dtype='object')
        temp_grids[:] = self.grids
        self.grids = temp_grids
        self.filenames = fn
        self._store_binary_hierarchy()
        t2 = time.time()

    def _fill_arrays(self, ei, si, LE, RE, np):
        self.grid_dimensions.flat[:] = ei
        self.grid_dimensions -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge.flat[:] = LE
        self.grid_right_edge.flat[:] = RE
        self.grid_particle_count.flat[:] = np

    def __pointer_handler(self, m):
        sgi = int(m[2])-1
        if sgi == -1: return # if it's 0, then we're done with that lineage
        # Okay, so, we have a pointer.  We make a new grid, with an id of the length+1
        # (recall, Enzo grids are 1-indexed)
        self.grids.append(self.grid(len(self.grids)+1, self))
        # We'll just go ahead and make a weakref to cache
        second_grid = self.grids[sgi] # zero-indexed already
        first_grid = self.grids[int(m[0])-1]
        if m[1] == "Next":
            first_grid._children_ids.append(second_grid.id)
            second_grid._parent_id = first_grid.id
            second_grid.Level = first_grid.Level + 1
        elif m[1] == "This":
            if first_grid.Parent is not None:
                first_grid.Parent._children_ids.append(second_grid.id)
                second_grid._parent_id = first_grid._parent_id
            second_grid.Level = first_grid.Level
        self.grid_levels[sgi] = second_grid.Level

    def _parse_binary_hierarchy(self):
        mylog.info("Getting the binary hierarchy")
        if not ytcfg.getboolean("yt","serialize"): return False
        try:
            f = h5py.File(self.hierarchy_filename[:-9] + "harrays")
        except:
            return False
        self.grid_dimensions[:] = f["/ActiveDimensions"][:]
        self.grid_left_edge[:] = f["/LeftEdges"][:]
        self.grid_right_edge[:] = f["/RightEdges"][:]
        self.grid_particle_count[:,0] = f["/NumberOfParticles"][:]
        levels = f["/Level"][:]
        parents = f["/ParentIDs"][:]
        procs = f["/Processor"][:]
        grids = []
        self.filenames = []
        grids = [self.grid(gi+1, self) for gi in xrange(self.num_grids)]
        giter = izip(grids, levels, procs, parents)
        bn = self._bn % (self.pf)
        pmap = [(bn % P,) for P in xrange(procs.max()+1)]
        for grid,L,P,Pid in giter:
            grid.Level = L
            grid._parent_id = Pid
            if Pid > -1:
                grids[Pid-1]._children_ids.append(grid.id)
            self.filenames.append(pmap[P])
        self.grids = na.array(grids, dtype='object')
        f.close()
        mylog.info("Finished with binary hierarchy reading")
        return True

    @parallel_blocking_call
    def _store_binary_hierarchy(self):
        # We don't do any of the logic here, we just check if the data file
        # is open...
        if self._data_file is None: return
        if self._data_mode == 'r': return
        if self.data_style != "enzo_packed_3d": return
        mylog.info("Storing the binary hierarchy")
        try:
            f = h5py.File(self.hierarchy_filename[:-9] + "harrays", "w")
        except IOError:
            return
        f.create_dataset("/LeftEdges", data=self.grid_left_edge)
        f.create_dataset("/RightEdges", data=self.grid_right_edge)
        parents, procs, levels = [], [], []
        for i,g in enumerate(self.grids):
            if g.Parent is not None:
                parents.append(g.Parent.id)
            else:
                parents.append(-1)
            procs.append(int(self.filenames[i][0][-4:]))
            levels.append(g.Level)

        parents = na.array(parents, dtype='int64')
        procs = na.array(procs, dtype='int64')
        levels = na.array(levels, dtype='int64')
        f.create_dataset("/ParentIDs", data=parents)
        f.create_dataset("/Processor", data=procs)
        f.create_dataset("/Level", data=levels)

        f.create_dataset("/ActiveDimensions", data=self.grid_dimensions)
        f.create_dataset("/NumberOfParticles", data=self.grid_particle_count[:,0])

        f.close()

    def _rebuild_top_grids(self, level = 0):
        #for level in xrange(self.max_level+1):
        mylog.info("Rebuilding grids on level %s", level)
        cmask = (self.grid_levels.flat == (level + 1))
        cmsum = cmask.sum()
        mask = na.zeros(self.num_grids, dtype='bool')
        for grid in self.select_grids(level):
            mask[:] = 0
            LE = self.grid_left_edge[grid.id - grid._id_offset]
            RE = self.grid_right_edge[grid.id - grid._id_offset]
            grids, grid_i = self.get_box_grids(LE, RE)
            mask[grid_i] = 1
            grid._children_ids = []
            cgrids = self.grids[ ( mask * cmask).astype('bool') ]
            mylog.info("%s: %s / %s", grid, len(cgrids), cmsum)
            for cgrid in cgrids:
                grid._children_ids.append(cgrid.id)
                cgrid._parent_id = grid.id
        mylog.info("Finished rebuilding")

    def _populate_grid_objects(self):
        for g,f in izip(self.grids, self.filenames):
            g._prepare_grid()
            g._setup_dx()
            g.set_filename(f[0])
            #if g.Parent is not None: g._guess_properties_from_parent()
        del self.filenames # No longer needed.
        self.max_level = self.grid_levels.max()

    def _detect_fields(self):
        self.field_list = []
        # Do this only on the root processor to save disk work.
        if self.comm.rank == 0 or self.comm.rank == None:
            field_list = self.get_data("/", "DataFields")
            if field_list is None:
                mylog.info("Gathering a field list (this may take a moment.)")
                field_list = set()
                random_sample = self._generate_random_grids()
                for grid in random_sample:
                    if not hasattr(grid, 'filename'): continue
                    try:
                        gf = self.io._read_field_names(grid)
                    except self.io._read_exception:
                        mylog.debug("Grid %s is a bit funky?", grid.id)
                        continue
                    mylog.debug("Grid %s has: %s", grid.id, gf)
                    field_list = field_list.union(gf)
        else:
            field_list = None
        field_list = self.comm.mpi_bcast(field_list)
        self.save_data(list(field_list),"/","DataFields",passthrough=True)
        self.field_list = list(field_list)

    def _generate_random_grids(self):
        if self.num_grids > 40:
            starter = na.random.randint(0, 20)
            random_sample = na.mgrid[starter:len(self.grids)-1:20j].astype("int32")
            # We also add in a bit to make sure that some of the grids have
            # particles
            gwp = self.grid_particle_count > 0
            if na.any(gwp) and not na.any(gwp[(random_sample,)]):
                # We just add one grid.  This is not terribly efficient.
                first_grid = na.where(gwp)[0][0]
                random_sample.resize((21,))
                random_sample[-1] = first_grid
                mylog.debug("Added additional grid %s", first_grid)
            mylog.debug("Checking grids: %s", random_sample.tolist())
        else:
            random_sample = na.mgrid[0:max(len(self.grids)-1,1)].astype("int32")
        return self.grids[(random_sample,)]

    def find_particles_by_type(self, ptype, max_num=None, additional_fields=None):
        """
        Returns a structure of arrays with all of the particles'
        positions, velocities, masses, types, IDs, and attributes for
        a particle type **ptype** for a maximum of **max_num**
        particles.  If non-default particle fields are used, provide
        them in **additional_fields**.
        """
        # Not sure whether this routine should be in the general HierarchyType.
        if self.grid_particle_count.sum() == 0:
            mylog.info("Data contains no particles.");
            return None
        if additional_fields is None:
            additional_fields = ['metallicity_fraction', 'creation_time',
                                 'dynamical_time']
        pfields = [f for f in self.field_list if f.startswith('particle_')]
        nattr = self.parameter_file['NumberOfParticleAttributes']
        if nattr > 0:
            pfields += additional_fields[:nattr]
        # Find where the particles reside and count them
        if max_num is None: max_num = 1e100
        total = 0
        pstore = []
        for level in range(self.max_level, -1, -1):
            for grid in self.select_grids(level):
                index = na.where(grid['particle_type'] == ptype)[0]
                total += len(index)
                pstore.append(index)
                if total >= max_num: break
            if total >= max_num: break
        result = None
        if total > 0:
            result = {}
            for p in pfields:
                result[p] = na.zeros(total, 'float64')
            # Now we retrieve data for each field
            ig = count = 0
            for level in range(self.max_level, -1, -1):
                for grid in self.select_grids(level):
                    nidx = len(pstore[ig])
                    if nidx > 0:
                        for p in pfields:
                            result[p][count:count+nidx] = grid[p][pstore[ig]]
                        count += nidx
                    ig += 1
                    if count >= total: break
                if count >= total: break
            # Crop data if retrieved more than max_num
            if count > max_num:
                for p in pfields:
                    result[p] = result[p][0:max_num]
        return result


class EnzoHierarchyInMemory(EnzoHierarchy):

    grid = EnzoGridInMemory
    _enzo = None

    @property
    def enzo(self):
        if self._enzo is None:
            import enzo
            self._enzo = enzo
        return self._enzo

    def __init__(self, pf, data_style = None):
        self.data_style = data_style
        self.float_type = 'float64'
        self.parameter_file = weakref.proxy(pf) # for _obtain_enzo
        self.float_type = self.enzo.hierarchy_information["GridLeftEdge"].dtype
        self.directory = os.getcwd()
        AMRHierarchy.__init__(self, pf, data_style)

    def _initialize_data_storage(self):
        pass

    def _count_grids(self):
        self.num_grids = self.enzo.hierarchy_information["GridDimensions"].shape[0]

    def _parse_hierarchy(self):
        self._copy_hierarchy_structure()
        mylog.debug("Copying reverse tree")
        reverse_tree = self.enzo.hierarchy_information["GridParentIDs"].ravel().tolist()
        # Initial setup:
        mylog.debug("Reconstructing parent-child relationships")
        grids = []
        # We enumerate, so it's 0-indexed id and 1-indexed pid
        self.filenames = ["-1"] * self.num_grids
        for id,pid in enumerate(reverse_tree):
            grids.append(self.grid(id+1, self))
            grids[-1].Level = self.grid_levels[id, 0]
            if pid > 0:
                grids[-1]._parent_id = pid
                grids[pid-1]._children_ids.append(grids[-1].id)
        self.max_level = self.grid_levels.max()
        mylog.debug("Preparing grids")
        self.grids = na.empty(len(grids), dtype='object')
        for i, grid in enumerate(grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = None
            grid._prepare_grid()
            grid.proc_num = self.grid_procs[i,0]
            self.grids[i] = grid
        mylog.debug("Prepared")

    def _initialize_grid_arrays(self):
        EnzoHierarchy._initialize_grid_arrays(self)
        self.grid_procs = na.zeros((self.num_grids,1),'int32')

    def _copy_hierarchy_structure(self):
        # Dimensions are important!
        self.grid_dimensions[:] = self.enzo.hierarchy_information["GridEndIndices"][:]
        self.grid_dimensions -= self.enzo.hierarchy_information["GridStartIndices"][:]
        self.grid_dimensions += 1
        self.grid_left_edge[:] = self.enzo.hierarchy_information["GridLeftEdge"][:]
        self.grid_right_edge[:] = self.enzo.hierarchy_information["GridRightEdge"][:]
        self.grid_levels[:] = self.enzo.hierarchy_information["GridLevels"][:]
        self.grid_procs = self.enzo.hierarchy_information["GridProcs"].copy()
        self.grid_particle_count[:] = self.enzo.hierarchy_information["GridNumberOfParticles"][:]

    def save_data(self, *args, **kwargs):
        pass

    _cached_field_list = None
    _cached_derived_field_list = None

    def _detect_fields(self):
        if self.__class__._cached_field_list is None:
            EnzoHierarchy._detect_fields(self)
            self.__class__._cached_field_list = self.field_list
        else:
            self.field_list = self.__class__._cached_field_list

    def _setup_derived_fields(self):
        if self.__class__._cached_derived_field_list is None:
            EnzoHierarchy._setup_derived_fields(self)
            self.__class__._cached_derived_field_list = self.derived_field_list
        else:
            self.derived_field_list = self.__class__._cached_derived_field_list

    def _generate_random_grids(self):
        my_rank = self.comm.rank
        my_grids = self.grids[self.grid_procs.ravel() == my_rank]
        if len(my_grids) > 40:
            starter = na.random.randint(0, 20)
            random_sample = na.mgrid[starter:len(my_grids)-1:20j].astype("int32")
            mylog.debug("Checking grids: %s", random_sample.tolist())
        else:
            random_sample = na.mgrid[0:max(len(my_grids)-1,1)].astype("int32")
        return my_grids[(random_sample,)]

class EnzoHierarchy1D(EnzoHierarchy):

    def _fill_arrays(self, ei, si, LE, RE, np):
        self.grid_dimensions[:,:1] = ei
        self.grid_dimensions[:,:1] -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge[:,:1] = LE
        self.grid_right_edge[:,:1] = RE
        self.grid_particle_count.flat[:] = np
        self.grid_left_edge[:,1:] = 0.0
        self.grid_right_edge[:,1:] = 1.0
        self.grid_dimensions[:,1:] = 1

class EnzoHierarchy2D(EnzoHierarchy):

    def _fill_arrays(self, ei, si, LE, RE, np):
        self.grid_dimensions[:,:2] = ei
        self.grid_dimensions[:,:2] -= na.array(si, self.float_type)
        self.grid_dimensions += 1
        self.grid_left_edge[:,:2] = LE
        self.grid_right_edge[:,:2] = RE
        self.grid_particle_count.flat[:] = np
        self.grid_left_edge[:,2] = 0.0
        self.grid_right_edge[:,2] = 1.0
        self.grid_dimensions[:,2] = 1

class EnzoStaticOutput(StaticOutput):
    """
    Enzo-specific output, set at a fixed time.
    """
    _hierarchy_class = EnzoHierarchy
    _fieldinfo_fallback = EnzoFieldInfo
    _fieldinfo_known = KnownEnzoFields
    def __init__(self, filename, data_style=None,
                 file_style = None,
                 parameter_override = None,
                 conversion_override = None,
                 storage_filename = None):
        """
        This class is a stripped down class that simply reads and parses
        *filename* without looking at the hierarchy.  *data_style* gets passed
        to the hierarchy to pre-determine the style of data-output.  However,
        it is not strictly necessary.  Optionally you may specify a
        *parameter_override* dictionary that will override anything in the
        paarmeter file and a *conversion_override* dictionary that consists
        of {fieldname : conversion_to_cgs} that will override the #DataCGS.
        """
        if filename.endswith(".hierarchy"): filename = filename[:-10]
        if parameter_override is None: parameter_override = {}
        self._parameter_override = parameter_override
        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override
        self.storage_filename = storage_filename

        StaticOutput.__init__(self, filename, data_style, file_style=file_style)
        if "InitialTime" not in self.parameters:
            self.current_time = 0.0
        rp = os.path.join(self.directory, "rates.out")
        if os.path.exists(rp):
            try:
                self.rates = EnzoTable(rp, rates_out_key)
            except:
                pass
        cp = os.path.join(self.directory, "cool_rates.out")
        if os.path.exists(cp):
            try:
                self.cool = EnzoTable(cp, cool_out_key)
            except:
                pass

        # Now fixes for different types of Hierarchies
        # This includes changing the fieldinfo class!
        if self["TopGridRank"] == 1: self._setup_1d()
        elif self["TopGridRank"] == 2: self._setup_2d()

    def _setup_1d(self):
        self._hierarchy_class = EnzoHierarchy1D
        self._fieldinfo_fallback = Enzo1DFieldInfo
        self.domain_left_edge = \
            na.concatenate([[self.domain_left_edge], [0.0, 0.0]])
        self.domain_right_edge = \
            na.concatenate([[self.domain_right_edge], [1.0, 1.0]])

    def _setup_2d(self):
        self._hierarchy_class = EnzoHierarchy2D
        self._fieldinfo_fallback = Enzo2DFieldInfo
        self.domain_left_edge = \
            na.concatenate([self["DomainLeftEdge"], [0.0]])
        self.domain_right_edge = \
            na.concatenate([self["DomainRightEdge"], [1.0]])

    def get_parameter(self,parameter,type=None):
        """
        Gets a parameter not in the parameterDict.
        """
        if self.parameters.has_key(parameter):
            return self.parameters[parameter]

        # Let's read the file
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        lines = open(self.parameter_filename).readlines()
        for lineI, line in enumerate(lines):
            if line.find("#") >= 1: # Keep the commented lines
                line=line[:line.find("#")]
            line=line.strip().rstrip()
            if len(line) < 2:
                continue
            try:
                param, vals = map(string.strip,map(string.rstrip,
                                                   line.split("=")))
            except ValueError:
                mylog.error("ValueError: '%s'", line)
            if parameter == param:
                if type is None:
                    t = vals.split()
                else:
                    t = map(type, vals.split())
                if len(t) == 1:
                    self.parameters[param] = t[0]
                else:
                    self.parameters[param] = t
                if param.endswith("Units") and not param.startswith("Temperature"):
                    dataType = param[:-5]
                    self.conversion_factors[dataType] = self.parameters[param]
                return self.parameters[parameter]

        return ""

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        # Let's read the file
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        lines = open(self.parameter_filename).readlines()
        data_labels = {}
        data_label_factors = {}
        for line in (l.strip() for l in lines):
            if len(line) < 2: continue
            param, vals = (i.strip() for i in line.split("="))
            # First we try to decipher what type of value it is.
            vals = vals.split()
            # Special case approaching.
            if "(do" in vals: vals = vals[:1]
            if len(vals) == 0:
                pcast = str # Assume NULL output
            else:
                v = vals[0]
                # Figure out if it's castable to floating point:
                try:
                    float(v)
                except ValueError:
                    pcast = str
                else:
                    if any("." in v or "e+" in v or "e-" in v for v in vals):
                        pcast = float
                    elif v == "inf":
                        pcast = str
                    else:
                        pcast = int
            # Now we figure out what to do with it.
            if param.endswith("Units") and not param.startswith("Temperature"):
                dataType = param[:-5]
                # This one better be a float.
                self.conversion_factors[dataType] = float(vals[0])
            if param.startswith("#DataCGS") or \
                 param.startswith("#CGSConversionFactor"):
                # Assume of the form: #DataCGSConversionFactor[7] = 2.38599e-26 g/cm^3
                # Which one does it belong to?
                data_id = param[param.find("[")+1:param.find("]")]
                data_label_factors[data_id] = float(vals[0])
            if param.startswith("DataLabel"):
                data_id = param[param.find("[")+1:param.find("]")]
                data_labels[data_id] = vals[0]
            if len(vals) == 0:
                vals = ""
            elif len(vals) == 1:
                vals = pcast(vals[0])
            else:
                vals = na.array([pcast(i) for i in vals if i != "-99999"])
            self.parameters[param] = vals
        for p, v in self._parameter_override.items():
            self.parameters[p] = v
        for p, v in self._conversion_override.items():
            self.conversion_factors[p] = v
        for k, v in data_label_factors.items():
            self.conversion_factors[data_labels[k]] = v
        self.refine_by = self.parameters["RefineBy"]
        self.dimensionality = self.parameters["TopGridRank"]
        if self.dimensionality > 1:
            self.domain_dimensions = self.parameters["TopGridDimensions"]
            if len(self.domain_dimensions) < 3:
                tmp = self.domain_dimensions.tolist()
                tmp.append(1)
                self.domain_dimensions = na.array(tmp)
            self.domain_left_edge = na.array(self.parameters["DomainLeftEdge"],
                                             "float64").copy()
            self.domain_right_edge = na.array(self.parameters["DomainRightEdge"],
                                             "float64").copy()
        else:
            self.domain_left_edge = na.array(self.parameters["DomainLeftEdge"],
                                             "float64")
            self.domain_right_edge = na.array(self.parameters["DomainRightEdge"],
                                             "float64")
            self.domain_dimensions = na.array([self.parameters["TopGridDimensions"],1,1])

        self.current_time = self.parameters["InitialTime"]
        # To be enabled when we can break old pickles:
        #if "MetaDataSimulationUUID" in self.parameters:
        #    self.unique_identifier = self.parameters["MetaDataSimulationUUID"]
        if "CurrentTimeIdentifier" in self.parameters:
            self.unique_identifier = self.parameters["CurrentTimeIdentifier"]
        if self.parameters["ComovingCoordinates"]:
            self.cosmological_simulation = 1
            self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
            self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
            self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
            self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        if "EOSType" not in self.parameters: self.parameters["EOSType"] = -1
        if self["ComovingCoordinates"]:
            self._setup_comoving_units()
        elif self.has_key("LengthUnit"):
            # 'Why share when we can reinvent incompatibly?'
            self.parameters["LengthUnits"] = self["LengthUnit"]
            self._setup_getunits_units()
        elif self.has_key("LengthUnits"):
            self._setup_getunits_units()
        else:
            self._setup_nounits_units()
        self.time_units['1'] = 1
        self.units['1'] = 1
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        self.time_units['Myr'] = self.time_units['years'] / 1.0e6
        self.time_units['Gyr']  = self.time_units['years'] / 1.0e9

    def _setup_comoving_units(self):
        z = self["CosmologyCurrentRedshift"]
        h = self["CosmologyHubbleConstantNow"]
        boxcm_cal = self["CosmologyComovingBoxSize"]
        boxcm_uncal = boxcm_cal / h
        box_proper = boxcm_uncal/(1+z)
        self.units['aye']  = (1.0 + self["CosmologyInitialRedshift"])/(z + 1.0)
        if not self.has_key("Time"):
            cu = self.cosmology_get_units()
            self.conversion_factors["Time"] = cu['utim']
        for unit in mpc_conversion:
            self.units[unit] = mpc_conversion[unit] * box_proper
            self.units[unit+'h'] = mpc_conversion[unit] * box_proper * h
            self.units[unit+'cm'] = mpc_conversion[unit] * boxcm_uncal
            self.units[unit+'hcm'] = mpc_conversion[unit] * boxcm_cal

    def _setup_getunits_units(self):
        # We are given LengthUnits, which is number of cm per box length
        # So we convert that to box-size in Mpc
        box_proper = 3.24077e-25 * self["LengthUnits"]
        self.units['aye']  = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] * box_proper
        if not self.has_key("TimeUnits"):
            self.conversion_factors["Time"] = self["LengthUnits"] / self["x-velocity"]

    def _setup_nounits_units(self):
        z = 0
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]

    def cosmology_get_units(self):
        """
        Return an Enzo-fortran style dictionary of units to feed into custom
        routines.  This is typically only necessary if you are interacting
        with fortran code.
        """
        k = {}
        k["utim"] = 2.52e17/na.sqrt(self.omega_matter)\
                       / self.hubble_constant \
                       / (1+self.parameters["CosmologyInitialRedshift"])**1.5
        k["urho"] = 1.88e-29 * self.omega_matter \
                        * self.hubble_constant**2 \
                        * (1.0 + self.current_redshift)**3
        k["uxyz"] = 3.086e24 * \
               self.parameters["CosmologyComovingBoxSize"] / \
               self.hubble_constant / \
               (1.0 + self.current_redshift)
        k["uaye"] = 1.0/(1.0 + self.parameters["CosmologyInitialRedshift"])
        k["uvel"] = 1.225e7*self.parameters["CosmologyComovingBoxSize"] \
                      *na.sqrt(self.omega_matter) \
                      *na.sqrt(1+ self.parameters["CosmologyInitialRedshift"])
        k["utem"] = 1.88e6 * (self.parameters["CosmologyComovingBoxSize"]**2) \
                      * self.omega_matter \
                      * (1.0 + self.parameters["CosmologyInitialRedshift"])
        k["aye"]  = (1.0 + self.parameters["CosmologyInitialRedshift"]) / \
               (1.0 + self.current_redshift)
        return k

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        if ("%s" % (args[0])).endswith(".hierarchy"):
            return True
        return os.path.exists("%s.hierarchy" % args[0])

class EnzoStaticOutputInMemory(EnzoStaticOutput):
    _hierarchy_class = EnzoHierarchyInMemory
    _data_style = 'enzo_inline'

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        obj.__init__(*args, **kwargs)
        return obj

    def __init__(self, parameter_override=None, conversion_override=None):
        if parameter_override is None: parameter_override = {}
        self._parameter_override = parameter_override
        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override

        StaticOutput.__init__(self, "InMemoryParameterFile", self._data_style)

    def _parse_parameter_file(self):
        enzo = self._obtain_enzo()
        self.basename = "cycle%08i" % (
            enzo.yt_parameter_file["NumberOfPythonCalls"])
        self.parameters['CurrentTimeIdentifier'] = time.time()
        self.parameters.update(enzo.yt_parameter_file)
        self.conversion_factors.update(enzo.conversion_factors)
        for i in self.parameters:
            if isinstance(self.parameters[i], types.TupleType):
                self.parameters[i] = na.array(self.parameters[i])
            if i.endswith("Units") and not i.startswith("Temperature"):
                dataType = i[:-5]
                self.conversion_factors[dataType] = self.parameters[i]
        self.domain_left_edge = self.parameters["DomainLeftEdge"].copy()
        self.domain_right_edge = self.parameters["DomainRightEdge"].copy()
        for i in self.conversion_factors:
            if isinstance(self.conversion_factors[i], types.TupleType):
                self.conversion_factors[i] = na.array(self.conversion_factors[i])
        for p, v in self._parameter_override.items():
            self.parameters[p] = v
        for p, v in self._conversion_override.items():
            self.conversion_factors[p] = v
        self.refine_by = self.parameters["RefineBy"]
        self.dimensionality = self.parameters["TopGridRank"]
        self.domain_dimensions = self.parameters["TopGridDimensions"]
        self.current_time = self.parameters["InitialTime"]
        if "CurrentTimeIdentifier" in self.parameters:
            self.unique_identifier = self.parameters["CurrentTimeIdentifier"]
        if self.parameters["ComovingCoordinates"]:
            self.cosmological_simulation = 1
            self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
            self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
            self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
            self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _obtain_enzo(self):
        import enzo; return enzo

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        return False

# These next two functions are taken from
# http://www.reddit.com/r/Python/comments/6hj75/reverse_file_iterator/c03vms4
# Credit goes to "Brian" on Reddit

def rblocks(f, blocksize=4096):
    """Read file as series of blocks from end of file to start.

    The data itself is in normal order, only the order of the blocks is reversed.
    ie. "hello world" -> ["ld","wor", "lo ", "hel"]
    Note that the file must be opened in binary mode.
    """
    if 'b' not in f.mode.lower():
        raise Exception("File must be opened using binary mode.")
    size = os.stat(f.name).st_size
    fullblocks, lastblock = divmod(size, blocksize)

    # The first(end of file) block will be short, since this leaves 
    # the rest aligned on a blocksize boundary.  This may be more 
    # efficient than having the last (first in file) block be short
    f.seek(-lastblock,2)
    yield f.read(lastblock)

    for i in range(fullblocks-1,-1, -1):
        f.seek(i * blocksize)
        yield f.read(blocksize)

def rlines(f, keepends=False):
    """Iterate through the lines of a file in reverse order.

    If keepends is true, line endings are kept as part of the line.
    """
    buf = ''
    for block in rblocks(f):
        buf = block + buf
        lines = buf.splitlines(keepends)
        # Return all lines except the first (since may be partial)
        if lines:
            lines.reverse()
            buf = lines.pop() # Last line becomes end of new first line.
            for line in lines:
                yield line
    yield buf  # First line.


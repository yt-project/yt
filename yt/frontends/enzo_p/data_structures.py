"""
Data structures for Enzo-P



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import \
    _h5py as h5py
import numpy as np
import os
import stat

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.static_output import \
    Dataset
from yt.fields.field_info_container import \
    NullFunc
from yt.funcs import \
    ensure_tuple, \
    get_pbar, \
    setdefaultattr
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.utilities.logger import \
    ytLogger as mylog

from yt.frontends.enzo_p.fields import \
    EnzoPFieldInfo
from yt.frontends.enzo_p.misc import \
    get_block_info, \
    get_child_index, \
    get_root_blocks, \
    get_root_block_id, \
    is_parent

class EnzoPGrid(AMRGridPatch):
    """
    Class representing a single EnzoP Grid instance.
    """

    _id_offset = 0
    _refine_by = 2

    def __init__(self, id, index, block_name, filename = None):
        """
        Returns an instance of EnzoPGrid with *id*, associated with
        *filename* and *index*.
        """
        # All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename=filename, index=index)
        self.block_name = block_name
        self._children_ids = None
        self._parent_id = -1
        self.Level = -1

    def __repr__(self):
        return "EnzoPGrid_%04d" % self.id

    def get_parent_id(self, desc_block_name):
        if self.block_name == desc_block_name:
            raise RuntimeError("Child and parent are the same!")
        dim = self.ds.dimensionality
        d_block = desc_block_name[1:].replace(":", "")
        parent = self

        while True:
            a_block = parent.block_name[1:].replace(":", "")
            gengap = (len(d_block) - len(a_block)) / dim
            if gengap <= 1:
                return parent.id
            cid = get_child_index(a_block, d_block)
            parent = self.index.grids[parent._children_ids[cid]]

    def add_child(self, child):
        if self._children_ids is None:
            self._children_ids = \
            -1*np.ones(self._refine_by**self.ds.dimensionality,
                       dtype=np.int64)

        a_block =  self.block_name[1:].replace(":", "")
        d_block = child.block_name[1:].replace(":", "")
        cid = get_child_index(a_block, d_block)
        self._children_ids[cid] = child.id

    @property
    def Parent(self):
        if self._parent_id == -1: return None
        return self.index.grids[self._parent_id]

    @property
    def Children(self):
        if self._children_ids is None:
            return []
        return [self.index.grids[cid]
                for cid in self._children_ids]

class EnzoPHierarchy(GridIndex):

    _strip_path = False
    grid = EnzoPGrid
    _preload_implemented = True

    def __init__(self, ds, dataset_type):

        self.dataset_type = dataset_type
        self.directory = os.path.dirname(ds.parameter_filename)
        self.index_filename = ds.parameter_filename
        if os.path.getsize(self.index_filename) == 0:
            raise IOError(-1,"File empty", self.index_filename)

        GridIndex.__init__(self, ds, dataset_type)
        self.dataset.dataset_type = self.dataset_type

    def _count_grids(self):
        fblock_size = 32768
        f = open(self.ds.parameter_filename, "r")
        f.seek(0, 2)
        file_size = f.tell()
        nblocks = np.ceil(float(file_size) /
                          fblock_size).astype(np.int64)
        f.seek(0)
        offset = f.tell()
        ngrids = 0
        for ib in range(nblocks):
            my_block = min(fblock_size, file_size - offset)
            buff = f.read(my_block)
            ngrids += buff.count("\n")
            offset += my_block
        f.close()
        self.num_grids = ngrids
        self.dataset_type = "enzo_p"

    def _parse_index(self):
        self.grids = np.empty(self.num_grids, dtype='object')

        pbar = get_pbar("Parsing Hierarchy ", self.num_grids)
        f = open(self.ds.parameter_filename, "r")
        fblock_size = 32768
        f.seek(0, 2)
        file_size = f.tell()
        nblocks = np.ceil(float(file_size) /
                          fblock_size).astype(np.int64)
        f.seek(0)
        offset = f.tell()
        lstr = ""
        # place child blocks after the root blocks
        rbdim = self.ds.root_block_dimensions
        nroot_blocks = rbdim.prod()
        child_id = nroot_blocks

        last_pid = None
        for ib in range(nblocks):
            fblock = min(fblock_size, file_size - offset)
            buff = lstr + f.read(fblock)
            bnl = 0
            for inl in range(buff.count("\n")):
                nnl = buff.find("\n", bnl)
                line = buff[bnl:nnl]
                block_name, block_file = line.split()
                level, left, right = get_block_info(block_name)

                rbindex = get_root_block_id(block_name)
                rbid = rbindex[0] * rbdim[1:].prod() + \
                  rbindex[1] * rbdim[2:].prod() + rbindex[2]
                if level == 0:
                    grid_id = rbid
                    parent_id = -1
                else:
                    grid_id = child_id
                    # Try the last parent_id first
                    if last_pid is not None and \
                      is_parent(self.grids[last_pid].block_name, block_name):
                        parent_id = last_pid
                    else:
                        parent_id = self.grids[rbid].get_parent_id(block_name)
                    last_pid = parent_id
                    child_id += 1

                my_grid = self.grid(
                    grid_id, self, block_name,
                    filename=os.path.join(self.directory, block_file))
                my_grid.Level = level
                my_grid._parent_id = parent_id

                self.grids[grid_id]           = my_grid
                self.grid_levels[grid_id]     = level
                self.grid_left_edge[grid_id]  = left
                self.grid_right_edge[grid_id] = right
                self.grid_dimensions[grid_id] = self.ds.active_grid_dimensions

                if level > 0:
                    self.grids[parent_id].add_child(my_grid)

                bnl = nnl + 1
                pbar.update(1)
            lstr = buff[bnl:]
            offset += fblock

        f.close()
        pbar.finish()

        slope = self.ds.domain_width / \
          self.ds.arr(np.ones(3), "code_length")
        self.grid_left_edge   = self.grid_left_edge  * slope + \
          self.ds.domain_left_edge
        self.grid_right_edge  = self.grid_right_edge * slope + \
          self.ds.domain_left_edge

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _setup_derived_fields(self):
        super(EnzoPHierarchy, self)._setup_derived_fields()
        for fname, field in self.ds.field_info.items():
            if not field.particle_type: continue
            if isinstance(fname, tuple): continue
            if field._function is NullFunc: continue

    def _detect_output_fields(self):
        self.field_list = []
        # Do this only on the root processor to save disk work.
        if self.comm.rank in (0, None):
            # Just check the first grid.
            grid = self.grids[0]
            field_list = self.io._read_field_names(grid)
            mylog.debug("Grid %s has: %s", grid.id, field_list)
            ptypes = self.dataset.particle_types
            ptypes_raw = self.dataset.particle_types_raw
        else:
            field_list = None
            ptypes = None
            ptypes_raw = None
        self.field_list = list(self.comm.mpi_bcast(field_list))
        self.dataset.particle_types = list(self.comm.mpi_bcast(ptypes))
        self.dataset.particle_types_raw = list(self.comm.mpi_bcast(ptypes_raw))

class EnzoPDataset(Dataset):
    """
    Enzo-P-specific output, set at a fixed time.
    """
    _index_class = EnzoPHierarchy
    _field_info_class = EnzoPFieldInfo

    def __init__(self, filename, dataset_type=None,
                 file_style = None,
                 parameter_override = None,
                 conversion_override = None,
                 storage_filename = None,
                 units_override=None,
                 unit_system="cgs"):
        """
        This class is a stripped down class that simply reads and parses
        *filename* without looking at the index.  *dataset_type* gets passed
        to the index to pre-determine the style of data-output.  However,
        it is not strictly necessary.  Optionally you may specify a
        *parameter_override* dictionary that will override anything in the
        paarmeter file and a *conversion_override* dictionary that consists
        of {fieldname : conversion_to_cgs} that will override the #DataCGS.
        """
        self.fluid_types += ("enzop",)
        if parameter_override is None: parameter_override = {}
        self._parameter_override = parameter_override
        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override
        self.storage_filename = storage_filename
        Dataset.__init__(self, filename, dataset_type, file_style=file_style,
                         units_override=units_override, unit_system=unit_system)

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """

        f = open(self.parameter_filename, "r")
        # get dimension from first block name
        b0, fn0 = f.readline().strip().split()
        level0, left0, right0 = get_block_info(b0, min_dim=0)
        root_blocks = get_root_blocks(b0)
        f.close()
        self.dimensionality = left0.size
        self.periodicity = \
          ensure_tuple(np.ones(self.dimensionality, dtype=bool))

        fh = h5py.File(os.path.join(self.directory, fn0), "r")
        self.domain_left_edge  = fh.attrs["lower"]
        self.domain_right_edge = fh.attrs["upper"]

        # all blocks are the same size
        ablock = fh[list(fh.keys())[0]]
        self.current_time = ablock.attrs["time"][0]
        gsi = ablock.attrs["enzo_GridStartIndex"]
        gei = ablock.attrs["enzo_GridEndIndex"]
        self.ghost_zones = gsi[0]
        self.root_block_dimensions = root_blocks
        self.active_grid_dimensions = gei - gsi + 1
        self.grid_dimensions = ablock.attrs["enzo_GridDimension"]
        self.domain_dimensions = root_blocks * self.active_grid_dimensions
        fh.close()

        self.periodicity += (False, ) * (3 - self.dimensionality)

        # WIP hard-coded for now
        self.refine_by = 2
        self.cosmological_simulation = 0
        self.gamma = 5. / 3.
        self.particle_types = ()
        self.particle_types_raw = self.particle_types
        self.unique_identifier = \
          str(int(os.stat(self.parameter_filename)[stat.ST_CTIME]))

    def _set_code_unit_attributes(self):
        setdefaultattr(self, 'length_unit', self.quan(1, "cm"))
        setdefaultattr(self, 'mass_unit', self.quan(1, "g"))
        setdefaultattr(self, 'time_unit', self.quan(1, "s"))
        setdefaultattr(self, 'velocity_unit',
                       self.length_unit / self.time_unit)
        magnetic_unit = np.sqrt(4*np.pi * self.mass_unit /
                                (self.time_unit**2 * self.length_unit))
        magnetic_unit = np.float64(magnetic_unit.in_cgs())
        setdefaultattr(self, 'magnetic_unit',
                       self.quan(magnetic_unit, "gauss"))

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        fn = args[0]
        ddir = os.path.dirname(fn)
        if not fn.endswith(".block_list"):
            return False
        try:
            with open(fn, "r") as f:
                block, block_file = f.readline().strip().split()
                get_block_info(block)
                if not os.path.exists(os.path.join(ddir, block_file)):
                    return False
        except Exception:
            return False
        return True

"""
GAMER-specific data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import stat
import numpy as np
import weakref

from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.file_handler import \
    HDF5FileHandler
from .fields import GAMERFieldInfo
from yt.testing import assert_equal



class GAMERGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id,
                              filename = index.index_filename,
                              index    = index)
        self.Parent   = None    # do NOT initialize Parent as []
        self.Children = []
        self.Level    = level

    def __repr__(self):
        return 'GAMERGrid_%09i (dimension = %s)' % (self.id, self.ActiveDimensions)


class GAMERHierarchy(GridIndex):
    grid                 = GAMERGrid
    _preload_implemented = True # since gamer defines "_read_chunk_data" in io.py
    
    def __init__(self, ds, dataset_type = 'gamer'):
        self.dataset_type     = dataset_type
        self.dataset          = weakref.proxy(ds)
        self.index_filename   = self.dataset.parameter_filename
        self.directory        = os.path.dirname(self.index_filename)
        self._handle          = ds._handle
        self.float_type       = 'float64' # fixed even when FLOAT8 is off
        self._particle_handle = ds._particle_handle
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # find all field names in the current dataset
        self.field_list = [ ('gamer', v) for v in self._handle['Data'].keys() ]
    
    def _count_grids(self):
        # count the total number of patches at all levels  
        self.num_grids = self.dataset.parameters['NPatch'].sum()
        
    def _parse_index(self):
        parameters       = self.dataset.parameters
        gid0             = 0
        grid_corner      = self._handle['Tree/Corner'].value
        convert2physical = self._handle['Tree/Corner'].attrs['Cvt2Phy']

        self.grid_dimensions    [:] = parameters['PatchSize']
        self.grid_particle_count[:] = 0

        for lv in range(0, parameters['NLevel']):
            num_grids_level = parameters['NPatch'][lv]
            if num_grids_level == 0: break

            patch_scale = parameters['PatchSize']*parameters['CellScale'][lv]

            # set the level and edge of each grid
            # (left/right_edge are YT arrays in code units)
            self.grid_levels.flat[ gid0:gid0 + num_grids_level ] = lv
            self.grid_left_edge[ gid0:gid0 + num_grids_level ] \
                = grid_corner[ gid0:gid0 + num_grids_level ]*convert2physical
            self.grid_right_edge[ gid0:gid0 + num_grids_level ] \
                = (grid_corner[ gid0:gid0 + num_grids_level ] + patch_scale)*convert2physical

            gid0 += num_grids_level

        # allocate all grid objects
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels.flat[i])

        # maximum level with patches (which can be lower than MAX_LEVEL)
        self.max_level = self.grid_levels.max()
        
    def _populate_grid_objects(self):
        son_list = self._handle["Tree/Son"].value

        for gid in range(self.num_grids):
            grid     = self.grids.flat[gid]
            son_gid0 = son_list[gid]

            # set up the parent-children relationship
            if son_gid0 >= 0:
                grid.Children = [ self.grids.flat[son_gid0+s] for s in range(8) ]

            for son_grid in grid.Children: son_grid.Parent = grid

            # set up other grid attributes
            grid._prepare_grid()
            grid._setup_dx()

        # validate the parent-children relationship in the debug mode
        if self.dataset._debug:
            self._validate_parent_children_relasionship()

    # for _debug mode only
    def _validate_parent_children_relasionship(self):
        mylog.info('Validating the parent-children relationship ...')

        father_list = self._handle["Tree/Father"].value

        for grid in self.grids:
            # parent->children == itself
            if grid.Parent is not None:
                assert grid.Parent.Children[0+grid.id%8] is grid, \
                       'Grid %d, Parent %d, Parent->Children %d' % \
                       (grid.id, grid.Parent.id, grid.Parent.Children[0].id)

            # children->parent == itself
            for c in grid.Children:
                assert c.Parent is grid, \
                       'Grid %d, Children %d, Children->Parent %d' % \
                       (grid.id, c.id, c.Parent.id)

            # all refinement grids should have parent 
            if grid.Level > 0:
                assert grid.Parent is not None and grid.Parent.id >= 0, \
                       'Grid %d, Level %d, Parent %d' % \
                       (grid.id, grid.Level, \
                        grid.Parent.id if grid.Parent is not None else -999)

            # parent index is consistent with the loaded dataset
            if grid.Level > 0:
                father_gid = father_list[grid.id]
                assert father_gid == grid.Parent.id, \
                       'Grid %d, Level %d, Parent_Found %d, Parent_Expect %d'%\
                       (grid.id, grid.Level, grid.Parent.id, father_gid)

            # edges between children and parent
            if len(grid.Children) > 0:
                assert_equal(grid.LeftEdge,  grid.Children[0].LeftEdge )
                assert_equal(grid.RightEdge, grid.Children[7].RightEdge)
        mylog.info('Check passed')
               

class GAMERDataset(Dataset):
    _index_class      = GAMERHierarchy
    _field_info_class = GAMERFieldInfo
    _handle           = None
    _debug            = False # debug mode for the GAMER frontend
    
    def __init__(self, filename,
                 dataset_type      = 'gamer',
                 storage_filename  = None,
                 particle_filename = None, 
                 units_override    = None,
                 unit_system       = "cgs"):

        if self._handle is not None: return

        self.fluid_types      += ('gamer',)
        self._handle           = HDF5FileHandler(filename)
        self.particle_filename = particle_filename

        if self.particle_filename is None:
            self._particle_handle = self._handle
        else:
            try:
                self._particle_handle = HDF5FileHandler(self.particle_filename)
            except:
                raise IOError(self.particle_filename)

        # currently GAMER only supports refinement by a factor of 2
        self.refine_by = 2

        Dataset.__init__(self, filename, dataset_type,
                         units_override = units_override,
                         unit_system    = unit_system)
        self.storage_filename = storage_filename
        
    def _set_code_unit_attributes(self):
        # GAMER does not assume any unit yet ...
        if len(self.units_override) == 0:
            mylog.warning("GAMER does not assume any unit ==> " +
                          "Use units_override to specify the units")

        for unit, cgs in [("length", "cm"), ("time", "s"), ("mass", "g")]:
            setdefaultattr(self, "%s_unit"%unit, self.quan(1.0, cgs))

            if len(self.units_override) == 0:
                mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
        
    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # shortcuts for different simulation information
        KeyInfo   = self._handle['Info']['KeyInfo']
        InputPara = self._handle['Info']['InputPara']
        Makefile  = self._handle['Info']['Makefile']
        SymConst  = self._handle['Info']['SymConst']

        # simulation time and domain
        self.current_time      = KeyInfo['Time'][0]
        self.dimensionality    = 3  # always 3D
        self.domain_left_edge  = np.array([0.,0.,0.], dtype='float64')
        self.domain_right_edge = KeyInfo['BoxSize'].astype('float64')
        self.domain_dimensions = KeyInfo['NX0'].astype('int64')

        # periodicity
        periodic         = InputPara['Opt__BC_Flu'][0] == 0
        self.periodicity = (periodic,periodic,periodic)

        # cosmological parameters
        if Makefile['Comoving']:
            self.cosmological_simulation = 1
            self.current_redshift        = 1.0/self.current_time - 1.0
            self.omega_matter            = InputPara['OmegaM0'] 
            self.omega_lambda            = 1.0 - self.omega_matter
            self.hubble_constant         = 0.6955   # H0 is not set in GAMER
        else:
            self.cosmological_simulation = 0
            self.current_redshift        = 0.0
            self.omega_matter            = 0.0
            self.omega_lambda            = 0.0
            self.hubble_constant         = 0.0

        # code-specific parameters
        for t in KeyInfo, InputPara, Makefile, SymConst:
            for v in t.dtype.names: self.parameters[v] = t[v]

        # reset 'Model' to be more readable
        if KeyInfo['Model'] == 1:
            self.parameters['Model'] = 'Hydro'
        elif KeyInfo['Model'] == 2:
            self.parameters['Model'] = 'MHD'
        elif KeyInfo['Model'] == 3:
            self.parameters['Model'] = 'ELBDM'
        else:
            self.parameters['Model'] = 'Unknown'

        # make aliases to some frequently used variables
        if self.parameters['Model'] == 'Hydro' or \
           self.parameters['Model'] == 'MHD':
            self.gamma = self.parameters["Gamma"]
            self.mu    = self.parameters.get("mu",0.6) # mean molecular weight

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            # define a unique way to identify GAMER datasets
            f = HDF5FileHandler(args[0])
            if 'Info' in f['/'].keys() and 'KeyInfo' in f['/Info'].keys():
                return True
        except:
            pass
        return False

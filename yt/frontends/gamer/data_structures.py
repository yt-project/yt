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

from yt.funcs import mylog
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
        self.PID      = None    # patch ID adopted in GAMER
        self.CID      = None    # cluster ID in the GAMER HDF5 output

    def __repr__(self):
        return 'GAMERGrid_%09i (dimension = %s)' % (self.id, self.ActiveDimensions)


class GAMERHierarchy(GridIndex):
    grid = GAMERGrid
    
    def __init__(self, ds, dataset_type = 'gamer'):
        self.dataset_type     = dataset_type
        self.dataset          = weakref.proxy(ds)
        self.index_filename   = self.dataset.parameter_filename
        self.directory        = os.path.dirname(self.index_filename)
        self._handle          = ds._handle
        self.float_type       = 'float64' # fixed even when FLOAT8 is off
#       self._particle_handle = ds._particle_handle
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # find all field names in the current dataset
        # choose an arbitrary grid since they all have the same fields
        Grid = self._handle['Level_00']['Cluster_%09i'%0]['Patch_%09i'%0]
        self.field_list = [ ('gamer', v) for v in Grid.dtype.names ]
    
    def _count_grids(self):
        # count the total number of patches at all levels  
        self.num_grids = self.dataset.parameters['NPatch'].sum()
        
    def _parse_index(self):
        f    = self._handle
        Para = self.dataset.parameters
        GID0 = 0

        self.grid_dimensions    [:] = Para['PatchSize']
        self.grid_particle_count[:] = 0

        for lv in range(0, Para['NLevel']):
            NP = Para['NPatch'][lv]
            if NP == 0: break

            Scale  = Para['CellScale'][lv]
            PScale = Para['PatchSize']*Scale
            CrList = f['PatchMap'][ 'Level_%d%d'%(lv/10,lv%10) ]['CornerList'].value
            Cr2Phy = f['PatchMap'][ 'Level_%d%d'%(lv/10,lv%10) ]['CornerList'].attrs['Cvt2Phy']

            # set the level and edge of each grid
            # (left/right_edge are YT arrays in code units)
            self.grid_levels.flat[ GID0:GID0+NP ] = lv
            self.grid_left_edge  [ GID0:GID0+NP ] = CrList[:]*Cr2Phy
            self.grid_right_edge [ GID0:GID0+NP ] = (CrList[:] + PScale)*Cr2Phy

            GID0 += NP

        # allocate all grid objects
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels.flat[i])

        # maximum level with patches (which can be lower than MAX_LEVEL)
        self.max_level = self.grid_levels.max()
        
    def _populate_grid_objects(self):
        f  = self._handle
        NP = self.dataset.parameters['NPatch']
        NC = f['SimuInfo']['KeyInfo']['H5_MaxDsetPerGroup']

        for lv in range(0, self.dataset.parameters['NLevel']):
            if NP[lv] == 0: break

            # converting HDF5 dataset to numpy array in advance (using .value)
            # is much faster than using h5py directly
            SonList    = f['PatchMap'][ 'Level_%02i'%lv ]['SonList'].value
#           SonList    = f['PatchMap'][ 'Level_%02i'%lv ]['SonList']
            SonGID0    = NP[0:lv+1].sum()
            GID0       = NP[0:lv  ].sum()

            # set the parent-children relationship
            for PID in range(0, NP[lv]):
                Grid   = self.grids.flat[GID0+PID]
                CID    = PID / NC
                SonPID = SonList[PID]

                if SonPID >= 0:
                    Grid.Children = [ self.grids.flat[SonGID0+SonPID+s] \
                                      for s in range(0,8) ]

                for SonGrid in Grid.Children: SonGrid.Parent = Grid

                # record the patch and cluster indices
                Grid.PID = PID
                Grid.CID = CID

                # set up other grid attributes
                Grid._prepare_grid()
                Grid._setup_dx()

        # validate the parent-children relationship in the debug mode
        if self.dataset._debug:
            self._validate_parent_children_relasionship()

    # for _debug mode only
    def _validate_parent_children_relasionship(self):
        mylog.info('Validating the parent-children relationship ...')

        f    = self._handle
        Para = self.dataset.parameters

        for Grid in self.grids:
            # parent->children == itself
            if Grid.Parent is not None:
                assert Grid.Parent.Children[0+Grid.id%8] is Grid, \
                       'Grid %d, Parent %d, Parent->Children %d' % \
                       (Grid.id, Grid.Parent.id, Grid.Parent.Children[0].id)

            # children->parent == itself
            for c in Grid.Children:
                assert c.Parent is Grid, \
                       'Grid %d, Children %d, Children->Parent %d' % \
                       (Grid.id, c.id, c.Parent.id)

            # all refinement grids should have parent 
            if Grid.Level > 0:
                assert Grid.Parent is not None and Grid.Parent.id >= 0, \
                       'Grid %d, Level %d, Parent %d' % \
                       (Grid.id, Grid.Level, \
                        Grid.Parent.id if Grid.Parent is not None else -999)

            # parent index is consistent with the loaded dataset
            if Grid.Level > 0:
                NC     = f['SimuInfo']['KeyInfo']['H5_MaxDsetPerGroup']
                PID    = Grid.PID
                CID    = PID / NC
                LvName = 'Level_%02i'   % Grid.Level
                CName  = 'Cluster_%09i' % CID
                PName  = 'Patch_%09i'   % PID
                FaGID  = f[LvName][CName][PName].attrs['Info']['Father'] \
                         + Para['NPatch'][0:Grid.Level-1].sum()
                assert FaGID == Grid.Parent.id, \
                       'Grid %d, Level %d, Parent_Found %d, Parent_Expect %d'%\
                       (Grid.id, Grid.Level, Grid.Parent.id, FaGID)

            # edges between children and parent
            if len(Grid.Children) > 0:
                assert_equal(Grid.LeftEdge,  Grid.Children[0].LeftEdge )
                assert_equal(Grid.RightEdge, Grid.Children[7].RightEdge)
        mylog.info('Check passed')
               

class GAMERDataset(Dataset):
    _index_class      = GAMERHierarchy
    _field_info_class = GAMERFieldInfo
    _handle           = None
#   _debug            = True  # turn on the debug mode
    _debug            = False # turn on the debug mode
    
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
            setattr(self, "%s_unit"%unit, self.quan(1.0, cgs))

            if len(self.units_override) == 0:
                mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
        
    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # shortcuts for different simulation information
        KeyInfo   = self._handle['SimuInfo']['KeyInfo']
        InputPara = self._handle['SimuInfo']['InputPara']
        Makefile  = self._handle['SimuInfo']['Makefile']
        SymConst  = self._handle['SimuInfo']['SymConst']

        # simulation time and domain
        self.current_time      = KeyInfo['Time'][0]
        self.dimensionality    = 3  # always 3D
        self.domain_left_edge  = np.array([0.,0.,0.], dtype='float64')
        self.domain_right_edge = KeyInfo['BoxSize'].astype('float64')
        self.domain_dimensions = KeyInfo['NX0'].astype('int64')

        # periodicity
        periodic = InputPara['Opt__BC_Flu'][0] == 0
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
        if KeyInfo['Model'] == 1:   self.parameters['Model'] = 'Hydro'
        elif KeyInfo['Model'] == 2: self.parameters['Model'] = 'MHD'
        elif KeyInfo['Model'] == 3: self.parameters['Model'] = 'ELBDM'
        else:                       self.parameters['Model'] = 'Unknown'

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
            if 'PatchMap' in f['/'].keys():
                return True
        except:
            pass
        return False

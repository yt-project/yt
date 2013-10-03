"""
TIGER-specific data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import *
from yt.data_objects.grid_patch import \
           AMRGridPatch
from yt.geometry.grid_geometry_handler import \
           GridIndex
from yt.data_objects.dataset import \
           Dataset

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from .fields import TigerFieldInfo, KnownTigerFields

class TigerGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, hierarchy, left_edge, right_edge, left_dims, right_dims):
        AMRGridPatch.__init__(self, id, hierarchy = hierarchy)
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        self.Level = 0
        self.NumberOfParticles = 0
        self.left_dims = np.array(left_dims, dtype='int32')
        self.right_dims = np.array(right_dims, dtype='int32')
        self.ActiveDimensions = self.right_dims - self.left_dims
        self.Parent = None
        self.Children = []

    @property
    def child_mask(self):
        return np.ones(self.ActiveDimensions, dtype='int32')

    def __repr__(self):
        return "TigerGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class TigerHierarchy(GridIndex):

    grid = TigerGrid

    def __init__(self, pf, dataset_type):
        self.directory = pf.fullpath
        self.dataset_type = dataset_type
        GridIndex.__init__(self, pf, dataset_type)

    def _count_grids(self):
        # Tiger is unigrid
        self.ngdims = [i/j for i,j in
                izip(self.pf.root_size, self.pf.max_grid_size)]
        self.num_grids = np.prod(self.ngdims)
        self.max_level = 0

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridIndex._setup_classes(self, dd)
        self.object_types.sort()

    def _parse_hierarchy(self):
        grids = []
        # We need to fill in dims, LE, RE, level, count
        dims, LE, RE, levels, counts = [], [], [], [], []
        DLE = self.pf.domain_left_edge
        DRE = self.pf.domain_right_edge 
        DW = DRE - DLE
        gds = DW / self.ngdims
        rd = [self.pf.root_size[i]-self.pf.max_grid_size[i] for i in range(3)]
        glx, gly, glz = np.mgrid[DLE[0]:DRE[0]-gds[0]:self.ngdims[0]*1j,
                                 DLE[1]:DRE[1]-gds[1]:self.ngdims[1]*1j,
                                 DLE[2]:DRE[2]-gds[2]:self.ngdims[2]*1j]
        gdx, gdy, gdz = np.mgrid[0:rd[0]:self.ngdims[0]*1j,
                                 0:rd[1]:self.ngdims[1]*1j,
                                 0:rd[2]:self.ngdims[2]*1j]
        LE, RE, levels, counts = [], [], [], []
        i = 0
        for glei, gldi in izip(izip(glx.flat, gly.flat, glz.flat),
                               izip(gdx.flat, gdy.flat, gdz.flat)):
            gld = np.array(gldi)
            gle = np.array(glei)
            gre = gle + gds
            g = self.grid(i, self, gle, gre, gld, gld+self.pf.max_grid_size)
            grids.append(g)
            dims.append(self.pf.max_grid_size)
            LE.append(g.LeftEdge)
            RE.append(g.RightEdge)
            levels.append(g.Level)
            counts.append(g.NumberOfParticles)
            i += 1
        self.grids = np.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self.grids[gi] = g
        self.grid_dimensions[:] = np.array(dims, dtype='int64')
        self.grid_left_edge[:] = np.array(LE, dtype='float64')
        self.grid_right_edge[:] = np.array(RE, dtype='float64')
        self.grid_levels.flat[:] = np.array(levels, dtype='int32')
        self.grid_particle_count.flat[:] = np.array(counts, dtype='int32')

    def _populate_grid_objects(self):
        # We don't need to do anything here
        for g in self.grids: g._setup_dx()

    def _detect_fields(self):
        self.file_mapping = {"Density" : "rhob",
                             "Temperature" : "temp"}

    @property
    def field_list(self):
        return self.file_mapping.keys()

    def _setup_derived_fields(self):
        self.derived_field_list = []

class TigerDataset(Dataset):
    _hierarchy_class = TigerHierarchy
    _fieldinfo_fallback = TigerFieldInfo
    _fieldinfo_known = KnownTigerFields

    def __init__(self, rhobname, root_size, max_grid_size=128,
                 dataset_type='tiger', storage_filename = None):
        Dataset.__init__(self, rhobname, dataset_type)
        self.storage_filename = storage_filename
        self.basename = rhobname[:-4]
        if not os.path.exists(self.basename + "rhob"):
            print "%s doesn't exist, don't know how to handle this!" % (
                        self.basename + "rhob")
            raise IOError
        if not iterable(root_size): root_size = (root_size,) * 3
        self.root_size = root_size
        if not iterable(max_grid_size): max_grid_size = (max_grid_size,) * 3
        self.max_grid_size = max_grid_size

        self.field_info = FieldInfoContainer.create_with_fallback(
                            self._fieldinfo_fallback)

        # We assume that we have basename + "rhob" and basename + "temp"
        # to get at our various parameters.

        # First we get our our header:
        
        header = [
            ('i', 'dummy0'),
            ('f', 'ZR'),
            ('f', 'OMEGA0'),
            ('f', 'FLAM0'),
            ('f', 'OMEGAB'),
            ('f', 'H0'),
            ('f', 'BOXL0'),
            ('i', 'dummy1'),
            ]

        h_fmt, h_key = zip(*header)
        header_string = "".join(h_fmt)

        fs = open(self.basename + "rhob")
        header_raw = read_struct(fs, header_string)
        self.parameters.update(dict(zip(h_key, header_raw)))

        if "InitialTime" not in self.parameters:
            self.current_time = 0.0
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[ST_CTIME])
        self.parameters['TopGridDimensions'] = root_size
        self.parameters['TopGridRank'] = 3
        self.units["Density"] = 1.0
        self.parameters['RefineBy'] = 2

    def _set_units(self):
        self.domain_left_edge = np.zeros(3, dtype='float64')
        self.domain_right_edge = np.ones(3, dtype='float64')
        self.units = {}
        self.time_units = {}
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['cm'] = 1.0 # This is just plain false
        self.units['unitary'] = 1.0 / (self["DomainRightEdge"] - self["DomainLeftEdge"]).max()

    def _parse_parameter_file(self):
        pass

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return os.path.exists(args[0] + "rhob")



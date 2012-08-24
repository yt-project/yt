"""
TIGER-specific data structures

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

from yt.funcs import *
from yt.data_objects.grid_patch import \
           AMRGridPatch
from yt.geometry.grid_geometry_handler import \
           GridGeometryHandler
from yt.data_objects.static_output import \
           StaticOutput

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

class TigerHierarchy(GridGeometryHandler):

    grid = TigerGrid

    def __init__(self, pf, data_style):
        self.directory = pf.fullpath
        self.data_style = data_style
        GridGeometryHandler.__init__(self, pf, data_style)

    def _count_grids(self):
        # Tiger is unigrid
        self.ngdims = [i/j for i,j in
                izip(self.pf.root_size, self.pf.max_grid_size)]
        self.num_grids = np.prod(self.ngdims)
        self.max_level = 0

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridGeometryHandler._setup_classes(self, dd)
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

class TigerStaticOutput(StaticOutput):
    _hierarchy_class = TigerHierarchy
    _fieldinfo_fallback = TigerFieldInfo
    _fieldinfo_known = KnownTigerFields

    def __init__(self, rhobname, root_size, max_grid_size=128,
                 data_style='tiger', storage_filename = None):
        StaticOutput.__init__(self, rhobname, data_style)
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



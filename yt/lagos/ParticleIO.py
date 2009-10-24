"""
The particle-IO handler

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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
from yt.lagos import *

particle_handler_registry = defaultdict()

class ParticleIOHandler(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_source_type"):
                particle_handler_registry[cls._source_type] = cls

    _source_type = None
    def __init__(self, pf, source):
        self.pf = pf
        self.data = {}
        self.source = source

    def __getitem__(self, key):
        if key not in self.data:
            self.get_data(key)
        return self.data[key]

    def __setitem__(self, key, val):
        self.data[key] = val

    def __delitem__(self, key):
        del self.data[key]

    def iter(self):
        for val in self.data.keys(): yield val

    def get_data(self, fields):
        fields = ensure_list(fields)
        self.source.get_data(fields)
        for field in fields:
            self[field] = self.source[field]

particle_handler_registry.default_factory = lambda: ParticleIOHandler

class ParticleIOHandlerRegion(ParticleIOHandler):
    periodic = False
    _source_type = "region"
    def __init__(self, pf, source):
        self.left_edge = source.left_edge
        self.right_edge = source.right_edge
        ParticleIOHandler.__init__(self, pf, source)

    def get_data(self, fields):
        mylog.info("Getting %s using ParticleIO" % str(fields))
        fields = ensure_list(fields)
        if not self.pf.h.io._particle_reader:
            mylog.info("not self.pf.h.io._particle_reader")
            return self.source.get_data(fields)
        rtype = 0
        DLE = na.array(self.pf["DomainLeftEdge"], dtype='float64') 
        DRE = na.array(self.pf["DomainRightEdge"], dtype='float64') 
        args = (na.array(self.left_edge), na.array(self.right_edge), 
                int(self.periodic), DLE, DRE)
        count_list, grid_list = [], []
        for grid in self.source._grids:
            if grid.NumberOfParticles == 0: continue
            grid_list.append(grid)
            if self.source._is_fully_enclosed(grid):
                count_list.append(grid.NumberOfParticles)
            else:
                count_list.append(-1)
        # region type, left_edge, right_edge, periodic, grid_list
        fields_to_read = []
        conv_factors = []
        for field in fields:
            f = self.pf.field_info[field]
            if f._particle_convert_function is None:
                fields_to_read.append(field)
                conv_factors.append(na.ones(len(grid_list), dtype='float64'))
            else:
                to_add = f.get_dependencies(pf = self.pf).requested
                if len(to_add) != 1: raise KeyError
                fields_to_read += to_add
                conv_factors.append(
                  na.fromiter((f.particle_convert(g) for g in grid_list),
                              count=len(grid_list), dtype='float64'))
        conv_factors = na.array(conv_factors).transpose()
        self.conv_factors = conv_factors
        rv = self.pf.h.io._read_particles(
            fields_to_read, rtype, args, grid_list, count_list,
            conv_factors)
        for field, v in zip(fields, rv): self[field] = v

class ParticleIOHandlerRegionStrict(ParticleIOHandlerRegion):
    _source_type = "region_strict"

class ParticleIOHandlerPeriodicRegion(ParticleIOHandlerRegion):
    periodic = True
    _source_type = "periodic_region"

class ParticleIOHandlerPeriodicRegionStrict(ParticleIOHandlerPeriodicRegion):
    _source_type = "periodic_region_strict"

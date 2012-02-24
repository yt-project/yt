"""
The particle-IO handler

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
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

import numpy as na

from yt.funcs import *

particle_handler_registry = defaultdict()

def particle_converter(func):
    from .data_containers import YTFieldData
    def save_state(grid):
        old_params = grid.field_parameters
        old_keys = grid.field_data.keys()
        tr = func(grid)
        grid.field_parameters = old_params
        grid.field_data = YTFieldData( [(k, grid.field_data[k]) for k in old_keys] )
        return tr
    return save_state

class ParticleIOHandler(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_source_type"):
                particle_handler_registry[cls._source_type] = cls

    _source_type = None

    def __init__(self, pf, source):
        self.pf = pf
        self.source = source

    def __getitem__(self, key):
        return self.get_data(key)

    def get_data(self, fields):
        fields = ensure_list(fields)
        rvs = self.source.get_data(fields, force_particle_read=True)
        if len(fields) == 1: return rvs[0]
        return rvs

particle_handler_registry.default_factory = lambda: ParticleIOHandler

class ParticleIOHandlerImplemented(ParticleIOHandler):
    def get_data(self, fields):
        mylog.info("Getting %s using ParticleIO" % str(fields))
        fields = ensure_list(fields)
        if not self.pf.h.io._particle_reader:
            mylog.info("not self.pf.h.io._particle_reader")
            return self.source.get_data(fields)
        rtype, args = self._get_args()
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
            to_add = f.get_dependencies(pf = self.pf).requested
            to_add = list(na.unique(to_add))
            if len(to_add) != 1: raise KeyError
            fields_to_read += to_add
            if f._particle_convert_function is None:
                func = f._convert_function
            else:
                func = f.particle_convert
            func = particle_converter(func)
            conv_factors.append(
              na.fromiter((func(g) for g in grid_list),
                          count=len(grid_list), dtype='float64'))
        conv_factors = na.array(conv_factors).transpose()
        self.conv_factors = conv_factors
        rvs = self.pf.h.io._read_particles(
            fields_to_read, rtype, args, grid_list, count_list,
            conv_factors)
        for [n, v] in zip(fields, rvs):
            self.source.field_data[n] = v

class ParticleIOHandlerRegion(ParticleIOHandlerImplemented):
    periodic = False
    _source_type = "region"

    def __init__(self, pf, source):
        self.left_edge = source.left_edge
        self.right_edge = source.right_edge
        ParticleIOHandler.__init__(self, pf, source)

    def _get_args(self):
        DLE = na.array(self.pf.domain_left_edge, dtype='float64') 
        DRE = na.array(self.pf.domain_right_edge, dtype='float64') 
        args = (na.array(self.left_edge), na.array(self.right_edge), 
                int(self.periodic), DLE, DRE)
        return (0, args)

class ParticleIOHandlerRegionStrict(ParticleIOHandlerRegion):
    _source_type = "region_strict"

class ParticleIOHandlerPeriodicRegion(ParticleIOHandlerRegion):
    periodic = True
    _source_type = "periodic_region"

class ParticleIOHandlerPeriodicRegionStrict(ParticleIOHandlerPeriodicRegion):
    _source_type = "periodic_region_strict"

class ParticleIOHandlerSphere(ParticleIOHandlerImplemented):
    _source_type = "sphere"

    def __init__(self, pf, source):
        self.center = source.center
        self.radius = source.radius
        ParticleIOHandler.__init__(self, pf, source)

    def _get_args(self):
        DLE = na.array(self.pf.domain_left_edge, dtype='float64')
        DRE = na.array(self.pf.domain_right_edge, dtype='float64')
        return (1, (na.array(self.center, dtype='float64'), self.radius,
            1, DLE, DRE))

class ParticleIOHandlerDisk(ParticleIOHandlerImplemented):
    _source_type = "disk"
    
    def __init__(self, pf, source):
        self.center = source.center
        self.normal = source._norm_vec
        self.radius = source._radius
        self.height = source._height
        ParticleIOHandler.__init__(self, pf, source)
    
    def _get_args(self):
        args = (na.array(self.center, dtype='float64'),
                na.array(self.normal, dtype='float64'),
                self.radius, self.height)
        return (2, args)
        

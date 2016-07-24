"""
The particle-IO handler



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from collections import defaultdict

from yt.funcs import \
    ensure_list, \
    mylog
from yt.extern.six import add_metaclass

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

class RegisteredParticleIOType(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if hasattr(cls, "_source_type"):
            particle_handler_registry[cls._source_type] = cls

@add_metaclass(RegisteredParticleIOType)
class ParticleIOHandler(object):
    _source_type = None

    def __init__(self, ds, source):
        self.ds = ds
        self.source = source

    def __getitem__(self, key):
        return self.get_data(key)

    def get_data(self, fields):
        fields = ensure_list(fields)
        self.source.get_data(fields, force_particle_read=True)
        rvs = [self.source[field] for field in fields]
        if len(fields) == 1: return rvs[0]
        return rvs

particle_handler_registry.default_factory = lambda: ParticleIOHandler

class ParticleIOHandlerImplemented(ParticleIOHandler):
    def get_data(self, fields):
        mylog.info("Getting %s using ParticleIO" % str(fields))
        fields = ensure_list(fields)
        if not self.ds.index.io._particle_reader:
            mylog.info("not self.ds.index.io._particle_reader")
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
            f = self.ds.field_info[field]
            to_add = f.get_dependencies(ds = self.ds).requested
            to_add = list(np.unique(to_add))
            if len(to_add) != 1: raise KeyError
            fields_to_read += to_add
            if f._particle_convert_function is None:
                func = f._convert_function
            else:
                func = f.particle_convert
            func = particle_converter(func)
            conv_factors.append(
              np.fromiter((func(g) for g in grid_list),
                          count=len(grid_list), dtype='float64'))
        conv_factors = np.array(conv_factors).transpose()
        self.conv_factors = conv_factors
        rvs = self.ds.index.io._read_particles(
            fields_to_read, rtype, args, grid_list, count_list,
            conv_factors)
        for [n, v] in zip(fields, rvs):
            self.source.field_data[n] = v

class ParticleIOHandlerRegion(ParticleIOHandlerImplemented):
    periodic = False
    _source_type = "region"

    def __init__(self, ds, source):
        self.left_edge = source.left_edge
        self.right_edge = source.right_edge
        ParticleIOHandler.__init__(self, ds, source)

    def _get_args(self):
        DLE = np.array(self.ds.domain_left_edge, dtype='float64') 
        DRE = np.array(self.ds.domain_right_edge, dtype='float64') 
        args = (np.array(self.left_edge), np.array(self.right_edge), 
                int(self.periodic), DLE, DRE)
        return (0, args)

class ParticleIOHandlerSphere(ParticleIOHandlerImplemented):
    _source_type = "sphere"

    def __init__(self, ds, source):
        self.center = source.center
        self.radius = source.radius
        ParticleIOHandler.__init__(self, ds, source)

    def _get_args(self):
        DLE = np.array(self.ds.domain_left_edge, dtype='float64')
        DRE = np.array(self.ds.domain_right_edge, dtype='float64')
        return (1, (np.array(self.center, dtype='float64'), self.radius,
            1, DLE, DRE))

class ParticleIOHandlerDisk(ParticleIOHandlerImplemented):
    _source_type = "disk"
    
    def __init__(self, ds, source):
        self.center = source.center
        self.normal = source._norm_vec
        self.radius = source.radius
        self.height = source.height
        ParticleIOHandler.__init__(self, ds, source)
    
    def _get_args(self):
        args = (np.array(self.center, dtype='float64'),
                np.array(self.normal, dtype='float64'),
                self.radius, self.height)
        return (2, args)
        

"""
The data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import h5py
import os
import re
import numpy as np
from yt.utilities.logger import ytLogger as mylog

from yt.utilities.io_handler import \
           BaseIOHandler

class IOHandlerChomboHDF5(BaseIOHandler):
    _dataset_type = "chombo_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, ds, *args, **kwargs):
        BaseIOHandler.__init__(self, ds, *args, **kwargs)
        self.ds = ds
        self._handle = ds._handle

    _field_dict = None
    @property
    def field_dict(self):
        if self._field_dict is not None:
            return self._field_dict
        field_dict = {}
        for key, val in self._handle.attrs.items():
            if key.startswith('component_'):
                comp_number = int(re.match('component_(\d)', key).groups()[0])
                field_dict[val] = comp_number
        self._field_dict = field_dict
        return self._field_dict

    _particle_field_index = None
    @property
    def particle_field_index(self):
        if self._particle_field_index is not None:
            return self._particle_field_index
        field_dict = {}
        for key, val in self._handle.attrs.items():
            if key.startswith('particle_'):
                comp_number = int(re.match('particle_component_(\d)', key).groups()[0])
                field_dict[val] = comp_number
        self._particle_field_index = field_dict
        return self._particle_field_index        
        
    def _read_field_names(self,grid):
        ncomp = int(self._handle.attrs['num_components'])
        fns = [c[1] for c in f.attrs.items()[-ncomp-1:-1]]
    
    def _read_data(self,grid,field):

        lstring = 'level_%i' % grid.Level
        lev = self._handle[lstring]
        dims = grid.ActiveDimensions
        boxsize = dims.prod()
        
        grid_offset = lev[self._offset_string][grid._level_id]
        start = grid_offset+self.field_dict[field]*boxsize
        stop = start + boxsize
        data = lev[self._data_string][start:stop]
        
        return data.reshape(dims, order='F')

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        chunks = list(chunks)
        fields.sort(key=lambda a: self.field_dict[a[1]])
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            grid = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[ftype, fname] = self._read_data(grid, fname)
            return rv
        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                   size, [f2 for f1, f2 in fields], ng)

        ind = 0
        for chunk in chunks:
            for g in chunk.objs:
                nd = 0
                for field in fields:
                    ftype, fname = field
                    data = self._read_data(g, fname)
                    nd = g.select(selector, data, rv[field], ind) # caches
                ind += nd
        return rv

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        chunks = list(chunks)

        if selector.__class__.__name__ == "GridSelector":

            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError

            grid = chunks[0].objs[0]

            for ftype, fname in fields:
                rv[ftype, fname] = self._read_particles(grid, fname)

            return rv

        rv = {f:np.array([]) for f in fields}
        for chunk in chunks:
            for grid in chunk.objs:
                for ftype, fname in fields:
                    data = self._read_particles(grid, fname)
                    rv[ftype, fname] = np.concatenate((data, rv[ftype, fname]))
        return rv

    def _read_particles(self, grid, name):

        field_index = self.particle_field_index[name]
        lev = 'level_%s' % grid.Level

        particles_per_grid = self._handle[lev]['particles:offsets'].value
        items_per_particle = len(self._particle_field_index)

        # compute global offset position
        offsets = items_per_particle * np.cumsum(particles_per_grid)
        offsets = np.append(np.array([0]), offsets)
        offsets = np.array(offsets, dtype=np.int64)

        # convert between the global grid id and the id on this level            
        grid_levels = np.array([g.Level for g in self.ds.index.grids])
        grid_ids    = np.array([g.id    for g in self.ds.index.grids])
        grid_level_offset = grid_ids[np.where(grid_levels == grid.Level)[0][0]]
        lo = grid.id - grid_level_offset
        hi = lo + 1

        # handle the case where this grid has no particles
        if (offsets[lo] == offsets[hi]):
            return np.array([], dtype=np.float64)

        data = self._handle[lev]['particles:data'][offsets[lo]:offsets[hi]]
        return np.asarray(data[field_index::items_per_particle], dtype=np.float64, order='F')

class IOHandlerChombo2DHDF5(IOHandlerChomboHDF5):
    _dataset_type = "chombo2d_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, ds, *args, **kwargs):
        BaseIOHandler.__init__(self, ds, *args, **kwargs)
        self.ds = ds
        self._handle = ds._handle

class IOHandlerChombo1DHDF5(IOHandlerChomboHDF5):
    _dataset_type = "chombo1d_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, ds, *args, **kwargs):
        BaseIOHandler.__init__(self, ds, *args, **kwargs)
        self.ds = ds
        self._handle = ds._handle   

class IOHandlerOrion2HDF5(IOHandlerChomboHDF5):
    _dataset_type = "orion_chombo_native"

    def _read_particles(self, grid, field):
        """
        parses the Orion Star Particle text files
             
        """

        fn = grid.ds.fullplotdir[:-4] + "sink"

        # Figure out the format of the particle file
        with open(fn, 'r') as f:
            lines = f.readlines()
        line = lines[1]

        # The basic fields that all sink particles have
        index = {'particle_mass': 0,
                 'particle_position_x': 1,
                 'particle_position_y': 2,
                 'particle_position_z': 3,
                 'particle_momentum_x': 4,
                 'particle_momentum_y': 5,
                 'particle_momentum_z': 6,
                 'particle_angmomen_x': 7,
                 'particle_angmomen_y': 8,
                 'particle_angmomen_z': 9,
                 'particle_id': -1}

        if len(line.strip().split()) == 11:
            # these are vanilla sinks, do nothing
            pass  

        elif len(line.strip().split()) == 17:
            # these are old-style stars, add stellar model parameters
            index['particle_mlast']     = 10
            index['particle_r']         = 11
            index['particle_mdeut']     = 12
            index['particle_n']         = 13
            index['particle_mdot']      = 14,
            index['particle_burnstate'] = 15

        elif len(line.strip().split()) == 18:
            # these are the newer style, add luminosity as well
            index['particle_mlast']     = 10
            index['particle_r']         = 11
            index['particle_mdeut']     = 12
            index['particle_n']         = 13
            index['particle_mdot']      = 14,
            index['particle_burnstate'] = 15,
            index['particle_luminosity']= 16

        else:
            # give a warning if none of the above apply:
            mylog.warning('Warning - could not figure out particle output file')
            mylog.warning('These results could be nonsense!')
            
        def read(line, field):
            return float(line.strip().split(' ')[index[field]])

        fn = grid.ds.fullplotdir[:-4] + "sink"
        with open(fn, 'r') as f:
            lines = f.readlines()
            particles = []
            for line in lines[1:]:
                if grid.NumberOfParticles > 0:
                    coord = read(line, "particle_position_x"), \
                            read(line, "particle_position_y"), \
                            read(line, "particle_position_z")
                    if ( (grid.LeftEdge < coord).all() and
                         (coord <= grid.RightEdge).all() ):
                        particles.append(read(line, field))
        return np.array(particles)

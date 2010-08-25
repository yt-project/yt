"""
The data-file handling functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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
from yt.lagos import *
import yt.amr_utils as au
import exceptions
import cPickle

_axis_ids = {0:2,1:1,2:0}

io_registry = {}

class BaseIOHandler(object):

    _data_style = None
    _particle_reader = False

    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_data_style"):
                io_registry[cls._data_style] = cls

    def __init__(self):
        self.queue = defaultdict(dict)

    # We need a function for reading a list of sets
    # and a function for *popping* from a queue all the appropriate sets

    def preload(self, grids, sets):
        pass

    def pop(self, grid, field):
        if grid.id in self.queue and field in self.queue[grid.id]:
            return self.modify(self.queue[grid.id].pop(field))
        else:
            # We only read the one set and do not store it if it isn't pre-loaded
            return self._read_data_set(grid, field)

    def peek(self, grid, field):
        return self.queue[grid.id].get(field, None)

    def push(self, grid, field, data):
        if grid.id in self.queue and field in self.queue[grid.id]:
            raise ValueError
        self.queue[grid][field] = data

    # Now we define our interface
    def _read_data_set(self, grid, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        pass

    def _read_field_names(self, grid):
        pass

    @property
    def _read_exception(self):
        return None

class IOHandlerExtracted(BaseIOHandler):

    _data_style = 'extracted'

    def _read_data_set(self, grid, field):
        return (grid.base_grid[field] / grid.base_grid.convert(field))

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return grid.base_grid[field][tuple(sl)] / grid.base_grid.convert(field)

class IOHandlerNative(BaseIOHandler):

    _data_style = "orion_native"

    def modify(self, field):
        return field.swapaxes(0,2)

    def _read_data_set(self,grid,field):
        """
        reads packed multiFABs output by BoxLib in "NATIVE" format.

        """
        filen = os.path.expanduser(grid.filename[field])
        off = grid._offset[field]
        inFile = open(filen,'rb')
        inFile.seek(off)
        header = inFile.readline()
        header.strip()

        if grid._paranoid:
            mylog.warn("Orion Native reader: Paranoid read mode.")
            headerRe = re.compile(orion_FAB_header_pattern)
            bytesPerReal,endian,start,stop,centerType,nComponents = headerRe.search(header).groups()

            # we will build up a dtype string, starting with endian
            # check endianness (this code is ugly. fix?)
            bytesPerReal = int(bytesPerReal)
            if bytesPerReal == int(endian[0]):
                dtype = '<'
            elif bytesPerReal == int(endian[-1]):
                dtype = '>'
            else:
                raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

            dtype += ('f%i'% bytesPerReal) #always a floating point

            # determine size of FAB
            start = na.array(map(int,start.split(',')))
            stop = na.array(map(int,stop.split(',')))

            gridSize = stop - start + 1

            error_count = 0
            if (start != grid.start).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid start." %grid.filename
                error_count += 1
            if (stop != grid.stop).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid stop." %grid.filename
                error_count += 1
            if (gridSize != grid.ActiveDimensions).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid dimensions." %grid.filename
                error_count += 1
            if bytesPerReal != grid.hierarchy._bytesPerReal:
                print "Paranoia Error: Cell_H and %s do not agree on bytes per real number." %grid.filename
                error_count += 1
            if (bytesPerReal == grid.hierarchy._bytesPerReal and dtype != grid.hierarchy._dtype):
                print "Paranoia Error: Cell_H and %s do not agree on endianness." %grid.filename
                error_count += 1

            if error_count > 0:
                raise RunTimeError("Paranoia unveiled %i differences between Cell_H and %s." % (error_count, grid.filename))

        else:
            start = grid.start_index
            stop = grid.stop_index
            dtype = grid.hierarchy._dtype
            bytesPerReal = grid.hierarchy._bytesPerReal

        nElements = grid.ActiveDimensions.prod()

        # one field has nElements*bytesPerReal bytes and is located
        # nElements*bytesPerReal*field_index from the offset location
        if yt2orionFieldsDict.has_key(field):
            fieldname = yt2orionFieldsDict[field]
        else:
            fieldname = field
        field_index = grid.field_indexes[fieldname]
        inFile.seek(int(nElements*bytesPerReal*field_index),1)
        field = na.fromfile(inFile,count=nElements,dtype=dtype)
        field = field.reshape(grid.ActiveDimensions[::-1]).swapaxes(0,2)

        # we can/should also check against the max and min in the header file

        inFile.close()
        return field

    def _read_data_slice(self, grid, field, axis, coord):
        """wishful thinking?
        """
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        #sl = tuple(reversed(sl))
        return self._read_data_set(grid,field)[sl]


class IOHandlerGadget(BaseIOHandler):
    _data_style = 'gadget_hdf5'
    def _read_data_set(self, grid, field):
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        if 'particles' in fh[adr].keys():
            adr2 = adr+'/particles'
            return fh[adr2][field]
        return None
    def _read_field_names(self,grid): 
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        rets = cPickle.loads(fh['/root'].attrs['fieldnames'])
        return rets

    def _read_data_slice(self,grid, field, axis, coord):
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        if 'particles' in fh[adr].keys():
            adr2 = adr+'/particles'
            return fh[adr2][field][coord,axis]
        return None
#
# Chombo readers
#

class IOHandlerChomboHDF5(BaseIOHandler):
    _data_style = "chombo_hdf5"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def _field_dict(self,fhandle):
        ncomp = int(fhandle['/'].attrs['num_components'])
        temp =  fhandle['/'].attrs.listitems()[-ncomp:]
        val, keys = zip(*temp)
        val = [int(re.match('component_(\d+)',v).groups()[0]) for v in val]
        return dict(zip(keys,val))
        
    def _read_field_names(self,grid):
        fhandle = h5py.File(grid.filename,'r')
        ncomp = int(fhandle['/'].attrs['num_components'])

        return [c[1] for c in f['/'].attrs.listitems()[-ncomp:]]
    
    def _read_data_set(self,grid,field):
        fhandle = h5py.File(grid.hierarchy.hierarchy_filename,'r')

        field_dict = self._field_dict(fhandle)
        lstring = 'level_%i' % grid.Level
        lev = fhandle[lstring]
        dims = grid.ActiveDimensions
        boxsize = dims.prod()
        
        grid_offset = lev[self._offset_string][grid._level_id]
        start = grid_offset+field_dict[field]*boxsize
        stop = start + boxsize
        data = lev[self._data_string][start:stop]

        return data.reshape(dims, order='F')
                                          

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid,field)[sl]

class IOHandlerTiger(BaseIOHandler):
    _data_style = "tiger"
    _offset = 36

    def __init__(self, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self._memmaps = {}

    def _read_data_set(self, grid, field):
        fn = grid.pf.basename + grid.hierarchy.file_mapping[field]
        LD = na.array(grid.left_dims, dtype='int64')
        SS = na.array(grid.ActiveDimensions, dtype='int64')
        RS = na.array(grid.pf.root_size, dtype='int64')
        data = au.read_tiger_section(fn, LD, SS, RS).astype("float64")
        return data

    def _read_data_slice(self, grid, field, axis, coord):
        fn = grid.pf.basename + grid.hierarchy.file_mapping[field]
        LD = na.array(grid.left_dims, dtype='int64').copy()
        SS = na.array(grid.ActiveDimensions, dtype='int64').copy()
        RS = na.array(grid.pf.root_size, dtype='int64').copy()
        LD[axis] += coord
        SS[axis] = 1
        data = au.read_tiger_section(fn, LD, SS, RS).astype("float64")
        return data

class IOHandlerFLASH(BaseIOHandler):
    _data_style = "flash_hdf5"

    def __init__(self, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)

    def _read_data_set(self, grid, field):
        f = h5py.File(grid.pf.parameter_filename, "r")
        tr = f["/%s" % field][grid.id - grid._id_offset,:,:,:].transpose()
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        f = h5py.File(grid.pf.parameter_filename, "r")
        tr = f["/%s" % field][grid.id - grid._id_offset].transpose()[sl]
        return tr.astype("float64")

class IOHandlerRAMSES(BaseIOHandler):
    _data_style = "ramses"

    def __init__(self, ramses_tree, *args, **kwargs):
        self.ramses_tree = ramses_tree
        BaseIOHandler.__init__(self, *args, **kwargs)

    def _read_data_set(self, grid, field):
        tr = na.zeros(grid.ActiveDimensions, dtype='float64')
        filled = na.zeros(grid.ActiveDimensions, dtype='int32')
        to_fill = grid.ActiveDimensions.prod()
        grids = [grid]
        l_delta = 0
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                to_fill -= self.ramses_tree.read_grid(field,
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, g.Level, 2**l_delta, g.locations)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid, field)[sl]

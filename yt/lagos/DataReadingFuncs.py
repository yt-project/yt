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
import exceptions

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

class IOHandlerHDF4(BaseIOHandler):

    _data_style = "enzo_hdf4"

    def modify(self, field):
        return field.swapaxes(0,2)

    def _read_field_names(self, grid):
        """
        Returns a list of fields associated with the filename
        Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
        """
        return SD.SD(grid.filename).datasets().keys()

    def _read_data_set(self, grid, field):
        """
        Returns after having obtained or generated a field.  Should throw an
        exception.  Should only be called as EnzoGridInstance.readData()

        @param field: field to read
        @type field: string
        """
        return SD.SD(grid.filename).select(field).get()

    def _read_data_slice(self, grid, field, axis, coord):
        """
        Reads a slice through the HDF4 data

        @param grid: Grid to slice
        @type grid: L{EnzoGrid<EnzoGrid>}
        @param field: field to get
        @type field: string
        @param sl: region to get
        @type sl: SliceType
        """
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        sl = tuple(reversed(sl))
        return SD.SD(grid.filename).select(field)[sl].swapaxes(0,2)

    @property
    def _read_exception(self):
        return SD.HDF4Error

class IOHandlerHDF4_2D(IOHandlerHDF4):

    _data_style = "enzo_hdf4_2d"

    def _read_data_set(self, grid, field):
        t = SD.SD(grid.filename).select(field).get()[:,:,None]
        return t.swapaxes(0,1)

    def _read_data_slice(self, grid, field, axis, coord):
        t = SD.SD(grid.filename).select(field).get()
        return t.transpose()

    def modify(self, field):
        return field

class IOHandlerHDF5(BaseIOHandler):

    _data_style = "enzo_hdf5"

    def _read_field_names(self, grid):
        """
        Returns a list of fields associated with the filename
        Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
        """
        return HDF5LightReader.ReadListOfDatasets(grid.filename, "/")

    def _read_data_set(self, grid, field):
        return HDF5LightReader.ReadData(grid.filename, "/%s" % field).swapaxes(0,2)

    def _read_data_slice(self, grid, field, axis, coord):
        """
        Reads a slice through the HDF5 data

        @param grid: Grid to slice
        @type grid: L{EnzoGrid<EnzoGrid>}
        @param field: field to get
        @type field: string
        @param axis: axis to slice along
        @param coord: coord to slice at
        """
        axis = {0:2,1:1,2:0}[axis]
        t = HDF5LightReader.ReadDataSlice(grid.filename, "/%s" %
                        (field), axis, coord).transpose()
        return t

    def modify(self, field):
        return field.swapaxes(0,2)

    @property
    def _read_exception(self):
        return (exceptions.KeyError, HDF5LightReader.ReadingError)


class IOHandlerPackedHDF5(BaseIOHandler):

    _data_style = "enzo_packed_3d"
    _particle_reader = True

    def _read_data_set(self, grid, field):
        return readDataPacked(grid, field)

    def _read_particles(self, fields, rtype, args, grid_list, enclosed):
        filenames = [g.filename for g in grid_list]
        ids = [g.id for g in grid_list]
        return HDF5LightReader.ReadParticles(
            rtype, fields, filenames, ids, args)

    def modify(self, field):
        return field.swapaxes(0,2)

    def preload(self, grids, sets):
        # We need to deal with files first
        files_keys = defaultdict(lambda: [])
        pf_field_list = grids[0].pf.h.field_list
        sets = [dset for dset in list(sets) if dset in pf_field_list]
        for g in grids: files_keys[g.filename].append(g)
        exc = self._read_exception
        for file in files_keys:
            mylog.debug("Starting read %s (%s)", file, sets)
            nodes = [g.id for g in files_keys[file]]
            nodes.sort()
            # We want to pass on any error we might expect -- the preload
            # phase should be non-fatal in all cases, and instead dump back to
            # the grids.
            data = HDF5LightReader.ReadMultipleGrids(file, nodes, sets)
            mylog.debug("Read %s items from %s", len(data), os.path.basename(file))
            for gid in data: self.queue[gid].update(data[gid])
        mylog.debug("Finished read of %s", sets)

    def _read_data_set(self, grid, field):
        return HDF5LightReader.ReadData(grid.filename,
                "/Grid%08i/%s" % (grid.id, field)).swapaxes(0,2)

    def _read_data_slice(self, grid, field, axis, coord):
        axis = _axis_ids[axis]
        return HDF5LightReader.ReadDataSlice(grid.filename, "/Grid%08i/%s" %
                        (grid.id, field), axis, coord).transpose()

    def _read_field_names(self, grid):
        return HDF5LightReader.ReadListOfDatasets(
                    grid.filename, "/Grid%08i" % grid.id)


class IOHandlerInMemory(BaseIOHandler):

    _data_style = "enzo_inline"

    def __init__(self, ghost_zones=3):
        import enzo
        self.enzo = enzo
        self.grids_in_memory = enzo.grid_data
        self.old_grids_in_memory = enzo.old_grid_data
        self.my_slice = (slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones))
        BaseIOHandler.__init__(self)

    def _read_data_set(self, grid, field):
        if grid.id not in self.grids_in_memory: raise KeyError
        return self.grids_in_memory[grid.id][field].swapaxes(0,2)[self.my_slice]
        # We don't do this, because we currently do not interpolate
        coef1 = max((grid.Time - t1)/(grid.Time - t2), 0.0)
        coef2 = 1.0 - coef1
        t1 = enzo.yt_parameter_file["InitialTime"]
        t2 = enzo.hierarchy_information["GridOldTimes"][grid.id]
        return (coef1*self.grids_in_memory[grid.id][field] + \
                coef2*self.old_grids_in_memory[grid.id][field])\
                [self.my_slice]

    def modify(self, field):
        return field.swapaxes(0,2)

    def _read_field_names(self, grid):
        return self.grids_in_memory[grid.id].keys()

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(3,-3), slice(3,-3), slice(3,-3)]
        sl[axis] = slice(coord + 3, coord + 4)
        sl = tuple(reversed(sl))
        return self.grids_in_memory[grid.id][field][sl].swapaxes(0,2)

    @property
    def _read_exception(self):
        return KeyError

class IOHandlerNative(BaseIOHandler):

    _data_style = "orion_native"

    def _read_data_set(self, grid, field):
        return readDataNative(grid, field)

    def modify(self, field):
        return field.swapaxes(0,2)

class IOHandlerPacked2D(IOHandlerPackedHDF5):

    _data_style = "enzo_packed_2d"

    def _read_data_set(self, grid, field):
        return HDF5LightReader.ReadData(grid.filename,
            "/Grid%08i/%s" % (grid.id, field)).transpose()[:,:,None]

    def modify(self, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        t = HDF5LightReader.ReadData(grid.filename, "/Grid%08i/%s" %
                        (grid.id, field)).transpose()
        return t


class IOHandlerPacked1D(IOHandlerPackedHDF5):

    _data_style = "enzo_packed_1d"

    def _read_data_set(self, grid, field):
        return HDF5LightReader.ReadData(grid.filename,
            "/Grid%08i/%s" % (grid.id, field)).transpose()[:,None,None]

    def modify(self, field):
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        t = HDF5LightReader.ReadData(grid.filename, "/Grid%08i/%s" %
                        (grid.id, field))
        return t

#
# BoxLib/Orion data readers follow
#
def readDataNative(self,field):
    """
    reads packed multiFABs output by BoxLib in "NATIVE" format.

    """
    filen = os.path.expanduser(self.filename[field])
    off = self._offset[field]
    inFile = open(filen,'rb')
    inFile.seek(off)
    header = inFile.readline()
    header.strip()

    if self._paranoid:
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
        if (start != self.start).any():
            print "Paranoia Error: Cell_H and %s do not agree on grid start." %self.filename
            error_count += 1
        if (stop != self.stop).any():
            print "Paranoia Error: Cell_H and %s do not agree on grid stop." %self.filename
            error_count += 1
        if (gridSize != self.ActiveDimensions).any():
            print "Paranoia Error: Cell_H and %s do not agree on grid dimensions." %self.filename
            error_count += 1
        if bytesPerReal != self.hierarchy._bytesPerReal:
            print "Paranoia Error: Cell_H and %s do not agree on bytes per real number." %self.filename
            error_count += 1
        if (bytesPerReal == self.hierarchy._bytesPerReal and dtype != self.hierarchy._dtype):
            print "Paranoia Error: Cell_H and %s do not agree on endianness." %self.filename
            error_count += 1

        if error_count > 0:
            raise RunTimeError("Paranoia unveiled %i differences between Cell_H and %s." % (error_count, self.filename))

    else:
        start = self.start
        stop = self.stop
        dtype = self.hierarchy._dtype
        bytesPerReal = self.hierarchy._bytesPerReal
        
    nElements = self.ActiveDimensions.prod()

    # one field has nElements*bytesPerReal bytes and is located
    # nElements*bytesPerReal*field_index from the offset location
    if yt2orionFieldsDict.has_key(field):
        fieldname = yt2orionFieldsDict[field]
    else:
        fieldname = field
    field_index = self.field_indexes[fieldname]
    inFile.seek(int(nElements*bytesPerReal*field_index),1)
    field = na.fromfile(inFile,count=nElements,dtype=dtype)
    field = field.reshape(self.ActiveDimensions[::-1]).swapaxes(0,2)

    # we can/should also check against the max and min in the header file
    
    inFile.close()
    return field
    
def readAllDataNative():
    pass

def readDataSliceNative(self, grid, field, axis, coord):
    """wishful thinking?
    """
    sl = [slice(None), slice(None), slice(None)]
    sl[axis] = slice(coord, coord + 1)
    #sl = tuple(reversed(sl))
    return grid.readDataFast(field)[sl]

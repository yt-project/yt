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

def getFieldsHDF4(self):
    """
    Returns a list of fields associated with the filename
    Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
    """
    return SD.SD(self.filename).datasets().keys()

def getFieldsHDF5(self):
    """
    Returns a list of fields associated with the filename
    Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
    """
    return HDF5LightReader.ReadListOfDatasets(self.filename, "/")

def readDataHDF4(self, field):
    """
    Returns after having obtained or generated a field.  Should throw an
    exception.  Should only be called as EnzoGridInstance.readData()

    @param field: field to read
    @type field: string
    """
    t = SD.SD(self.filename).select(field).get()
    return t.swapaxes(0,2)

def readAllDataHDF4(self):
    """
    Reads all fields inside an HDF4 file.  Should only be called as
    EnzoGridInstance.readAllData() .
    """
    sets = SD.SD(self.filename).datasets()
    for set in sets:
        self[set] = self.readDataFast(set)

def readDataHDF5(self, field):
    return HDF5LightReader.ReadData(self.filename, "/%s" % field).swapaxes(0,2)

def readAllDataHDF5(self):
    """
    Not implemented.  Fix me!
    """
    pass

def readAllDataPacked(self):
    """
    Not implemented.  Fix me!
    """
    pass

def readDataSliceHDF5(self, grid, field, axis, coord):
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

def readDataSliceHDF4(self, grid, field, axis, coord):
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

def readDataPackedHandle(self, field):
    t = self.handle.getNode("/Grid%08i" % (self.id), field).read().astype('float64')
    t = t.swapaxes(0,2)
    return t

def readDataPacked(self, field):
    return HDF5LightReader.ReadData(self.filename, "/Grid%08i/%s" % (self.id, field)).swapaxes(0,2)

def readDataSlicePacked(self, grid, field, axis, coord):
    """
    Reads a slice through the HDF5 data

    @param grid: Grid to slice
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param field: field to get
    @type field: string
    @param sl: region to get
    @type sl: SliceType
    """
    axis = {0:2,1:1,2:0}[axis]
    t = HDF5LightReader.ReadDataSlice(grid.filename, "/Grid%08i/%s" %
                    (grid.id, field), axis, coord).transpose()
    return t

def getFieldsPacked(self):
    """
    Returns a list of fields associated with the filename
    Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
    """
    return HDF5LightReader.ReadListOfDatasets(self.filename, "/Grid%08i" % self.id)

def getExceptionHDF4():
    return SD.HDF4Error

def getExceptionHDF5():
    return (exceptions.KeyError, HDF5LightReader.ReadingError)

def readDataInMemory(self, field):
    import enzo
    return enzo.grid_data[self.id][field].swapaxes(0,2)

def readAllDataInMemory(self):
    pass

def getFieldsInMemory(self):
    import enzo
    return enzo.grid_data[self.id].keys()

def readDataSliceInMemory(self, grid, field, axis, coord):
    import enzo
    sl = [slice(3,-3), slice(3,-3), slice(3,-3)]
    sl[axis] = slice(coord + 3, coord + 4)
    sl = tuple(reversed(sl))
    return enzo.grid_data[grid.id][field][sl].swapaxes(0,2)

def getExceptionInMemory():
    return KeyError

def readDataSliceHDF4_2D(self, grid, field, axis, coord):
    t = SD.SD(grid.filename).select(field).get()
    return t.transpose()


def readDataSlicePacked2D(self, grid, field, axis, coord):
    """
    Reads a slice through the HDF5 data

    @param grid: Grid to slice
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param field: field to get
    @type field: string
    @param sl: region to get
    @type sl: SliceType
    """
    t = HDF5LightReader.ReadData(grid.filename, "/Grid%08i/%s" %
                    (grid.id, field)).transpose()
    return t

def readDataSlicePacked1D(self, grid, field, axis, coord):
    """
    Reads a slice through the HDF5 data

    @param grid: Grid to slice
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param field: field to get
    @type field: string
    @param sl: region to get
    @type sl: SliceType
    """
    t = HDF5LightReader.ReadData(grid.filename, "/Grid%08i/%s" %
                    (grid.id, field))
    return t

class BaseDataQueue(object):

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
            return self._read_set(grid, field)

    def peek(self, grid, field):
        return self.queue[grid.id].get(field, None)

    def push(self, grid, field, data):
        if grid.id in self.queue and field in self.queue[grid.id]:
            raise ValueError
        self.queue[grid][field] = data

class DataQueueHDF4(BaseDataQueue):
    def _read_set(self, grid, field):
        return readDataHDF4(grid, field)

    def modify(self, field):
        return field.swapaxes(0,2)

class DataQueueHDF4_2D(BaseDataQueue):
    def _read_set(self, grid, field):
        t = SD.SD(grid.filename).select(field).get()[:,:,None]
        return t.swapaxes(0,1)

    def modify(self, field):
        pass

class DataQueueHDF5(BaseDataQueue):
    def _read_set(self, grid, field):
        return readDataHDF5(grid, field)

    def modify(self, field):
        return field.swapaxes(0,2)

class DataQueuePackedHDF5(BaseDataQueue):
    def _read_set(self, grid, field):
        return readDataPacked(grid, field)

    def modify(self, field):
        return field.swapaxes(0,2)

    def preload(self, grids, sets):
        # We need to deal with files first
        files_keys = defaultdict(lambda: [])
        pf_field_list = grids[0].pf.h.field_list
        sets = [dset for dset in list(sets) if dset in pf_field_list]
        for g in grids: files_keys[g.filename].append(g)
        exc = getExceptionHDF5()
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

class DataQueueInMemory(BaseDataQueue):
    def __init__(self, ghost_zones=3):
        import enzo
        self.grids_in_memory = enzo.grid_data
        self.old_grids_in_memory = enzo.old_grid_data
        self.my_slice = (slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones),
                      slice(ghost_zones,-ghost_zones))
        BaseDataQueue.__init__(self)

    def _read_set(self, grid, field):
        import enzo
        if grid.id not in self.grids_in_memory: raise KeyError
        return self.grids_in_memory[grid.id][field][self.my_slice]
        coef1 = max((grid.Time - t1)/(grid.Time - t2), 0.0)
        coef2 = 1.0 - coef1
        t1 = enzo.yt_parameter_file["InitialTime"]
        t2 = enzo.hierarchy_information["GridOldTimes"][grid.id]
        return (coef1*self.grids_in_memory[grid.id][field] + \
                coef2*self.old_grids_in_memory[grid.id][field])\
                [self.my_slice]

    def modify(self, field):
        return field.swapaxes(0,2)

    def preload(self, grids, sets):
        pass

class DataQueueNative(BaseDataQueue):
    def _read_set(self, grid, field):
        return readDataNative(grid, field)

    def modify(self, field):
        return field.swapaxes(0,2)

class DataQueuePacked2D(BaseDataQueue):
    def _read_set(self, grid, field):
        return HDF5LightReader.ReadData(grid.filename,
            "/Grid%08i/%s" % (grid.id, field)).transpose()[:,:,None]

    def modify(self, field):
        pass

class DataQueuePacked1D(BaseDataQueue):
    def _read_set(self, grid, field):
        return HDF5LightReader.ReadData(grid.filename,
            "/Grid%08i/%s" % (grid.id, field)).transpose()[:,None,None]

    def modify(self, field):
        pass

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

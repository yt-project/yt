"""
The data-file handling functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

class BaseDataQueue(object):

    def __init__(self):
        self.queue = defaultdict(lambda: {})

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

class DataQueuePackedHDF5(BaseDataQueue):

    def _read_set(self, grid, field):
        return readDataPacked(grid, field)

    def modify(self, field):
        return field.swapaxes(0,2)

    def preload(self, grids, sets):
        # We need to deal with files first
        files_keys = defaultdict(lambda: [])
        sets = list(sets)
        for g in grids: files_keys[g.filename].append(g)
        for file in files_keys:
            mylog.debug("Starting read %s", file)
            nodes = [g.id for g in files_keys[file]]
            nodes.sort()
            data = HDF5LightReader.ReadMultipleGrids(file, nodes, sets)
            mylog.debug("Read %s items from %s", len(data), os.path.basename(file))
            for gid in data: self.queue[gid].update(data[gid])
        mylog.debug("Finished read of %s", sets)

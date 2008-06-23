"""
The data-file handling functions

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

def readDataNative(self,field):
    """
    reads packed multiFABs output by BoxLib in "NATIVE" format.

    """
    inFile = open(os.path.expanduser(self.filename),'rb')
    inFile.seek(self._offset)
    header = inFile.readline()
    header.strip()

    if self._paranoid:
        mylog.warn("Orion Native reader: Paranoid read mode.")
        pattern = r"^FAB \(\((\d+), \([0-9 ]+\)\),\(\d+, \(([0-9 ]+)\)\)\)\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n"

        headerRe = re.compile(pattern)
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

        self._dtype = dtype
        # determine size of FAB
        # TODO: we should check consistency of this against the MF header...
        start = na.array(map(int,start.split(',')))
        stop = na.array(map(int,stop.split(',')))

        gridSize = stop - start + 1
        if (gridSize != self.ActiveDimensions).any():
            pass
            #raise KeyError("Your paranoia was well warrented. Cell_H and %s do not agree on grid size." % self.filename)
    nElements = self.ActiveDimensions.prod()

    # one field has nElements*bytesPerReal bytes and is located
    # nElements*bytesPerReal*field_index from the offset location
    field_index = self.field_indexes[yt2orionFieldsDict[field]]
    inFile.seek(int(nElements*bytesPerReal*field_index),1)
    field = na.fromfile(inFile,count=nElements,dtype=self._dtype)
    field = field.reshape(self.ActiveDimensions[::-1]).swapaxes(0,2)

    # we can/should also check against the max and min in the header file
    
    inFile.close()
    return field
    return na.ones(self.ActiveDimensions, dtype='float64')#field

def readAllDataNative():
    pass

def readDataSliceNative(self, grid, field, axis, coord):
    """wishful thinking?
    """
    sl = [slice(None), slice(None), slice(None)]
    sl[axis] = slice(coord, coord + 1)
    sl = tuple(reversed(sl))
    return self.readDataFast(field)[sl]

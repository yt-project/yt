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
    fls = []
    file = tables.openFile(self.filename)
    for fl in file.listNodes("/"):
        fls.append(fl.name)
    file.close()
    return fls

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
        self.readDataFast(set)

def readDataHDF5(self, field):
    """
    Reads a field from an HDF5 file.  Should only be called as
    EnzoGridInstance.readData()

    @param field: field to read
    @type field: string
    """
    f = tables.openFile(self.filename)
    t = f.getNode("/", field).read().astype("float64")
    t = t.swapaxes(0,2)
    f.close()
    return t

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

def readDataSliceHDF5(self, grid, field, sl):
    """
    Reads a slice through the HDF5 data

    @param grid: Grid to slice
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param field: field to get
    @type field: string
    @param sl: region to get
    @type sl: SliceType
    """
    f = tables.openFile(grid.filename)
    ss = f.getNode("/", field)[sl].swapaxes(0,2)
    f.close()
    return ss

def readDataSliceHDF4(self, grid, field, sl):
    """
    Reads a slice through the HDF4 data

    @param grid: Grid to slice
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param field: field to get
    @type field: string
    @param sl: region to get
    @type sl: SliceType
    """
    return SD.SD(grid.filename).select(field)[sl].swapaxes(0,2)

def readDataPacked(self, field):
    f = tables.openFile(self.filename,
                        rootUEP="/Grid%08i" % (self.id),
                        mode='r')
    t = f.getNode("/", field).read()
    t = t.swapaxes(0,2)
    f.close()
    return t

def readDataSlicePacked(self, grid, field, sl):
    """
    Reads a slice through the HDF5 data

    @param grid: Grid to slice
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param field: field to get
    @type field: string
    @param sl: region to get
    @type sl: SliceType
    """
    f = tables.openFile(grid.filename)
    ss = f.getNode("/Grid%08i" % (grid.id), field)[sl]
    f.close()
    return ss

def getFieldsPacked(self):
    """
    Returns a list of fields associated with the filename
    Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
    """
    fls = []
    f = tables.openFile(self.filename,
                        rootUEP="/Grid%08i" % (self.id),
                        mode='r')
    for fl in f.listNodes("/"):
        fls.append(fl.name)
    f.close()
    del f
    return fls

def getExceptionHDF4():
    return SD.HDF4Error

def getExceptionHDF5():
    return exceptions.KeyError

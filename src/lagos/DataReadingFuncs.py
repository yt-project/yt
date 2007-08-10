"""
The data-file handling functions

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""
from yt.lagos import *

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
    for fl in tables.openFile(self.filename).listNodes("/"):
        fls.append(fl.name)
    return fls

def readDataHDF4(self, field):
    """
    Returns after having obtained or generated a field.  Should throw an
    exception.  Should only be called as EnzoGridInstance.readData()
    
    @param field: field to read
    @type field: string
    """
    if self.data.has_key(field):
        return 1
    try:
        t = SD.SD(self.filename).select(field).get()
        self[field] = t.swapaxes(0,2)
    except:
        self.generateField(field)
    return 2

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
    EnzoGridInstance.realData()

    @param field: field to read
    @type field: string
    """
    if self.has_key(field):
        return 1
    f = tables.openFile(self.filename)
    try:
        t = f.getNode("/", field).read().astype("float64")
        self[field] = t.swapaxes(0,2)
    except:
        self.generateField(field)
    #self[field] = ones(self.data[field].shape)
    f.close()
    return 2

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
    ss = f.getNode("/", field)[sl]
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
    return SD.SD(grid.filename).select(field)[sl]

def readDataPacked(self, field):
    if self.has_key(field):
        return 1
    f = tables.openFile(self.filename)
    try:
        t = f.getNode("/Grid%08i" % (self.id), field).read()
        self[field] = t.swapaxes(0,2)
    except:
        self.generateField(field)
    #self[field] = ones(self.data[field].shape)
    f.close()
    return 2

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
    f = tables.openFile(self.filename)
    for fl in f.listNodes("/Grid%08i" % (self.id)):
        fls.append(fl.name)
    f.close()
    return fls


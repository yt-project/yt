from yt.lagos import *

def getFieldsHDF4(self):
    """
    Returns a list of fields associated with the filename
    Should *only* be called as EnzoGridInstance.getFields, never as getFields(object)
    """
    return SD.SD(self.filename).datasets()

def getFieldsHDF5(self):
    """
    Current unimplemented.  Fix me, someone!
    """
    return 0

def readDataHDF4(self, field):
    """
    Returns after having obtained or generated a field.  Should throw an
    exception.  Should only be called as EnzoGridInstance.readData()
    
    Arguments:
        field -- field to read
    """
    if self.data.has_key(field):
        return 1
    try:
        self[field] = SD.SD(self.filename).select(field).get()
        self[field].swapaxes(0,2)
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

    Arguments:
        field -- field to read
    """
    if self.has_key(field):
        return 1
    f = tables.openFile(self.filename)
    try:
        self[field] = f.getNode("/", field).read()
        self[field].swapaxes(0,2)
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

def readDataSliceHDF4(self, grid, field, sl):
    return SD.SD(grid.filename).select(field)[sl]

"""
Module for writing out Enzo data.  May never be finished.

Basically, now that we have wrapped WriteAllData, this module is not
necessarily necessary.  It might be worthwhile to do spot-changes, though, of
existing grids.

@deprecated: Use the EnzoInterface module

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from pyhdf import SD
import pyhdf.error
import tables, warnings
from numarray import *
import numarray.objects as obj
import numarray.nd_image as nd
import numarray
from string import strip, rstrip
from math import ceil, floor, log10, pi
import os.path, types, exceptions

from ravenDefs import *

import fields

def writeFieldHDF4(grid, fieldName, newField):
    if fieldName in grid.data.keys():
        del grid.data[fieldName]
    s = SD.SD(grid.fileName)
    p = s.get(fieldName)
    p[:] = newField
    p.endaccess()
    s.end()



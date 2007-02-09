#
# write:
#   A module for modifying and writing out Enzo data
#
# Written by: Matthew Turk (mturk@stanford.edu) Dec 2006
# Modified:
#

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



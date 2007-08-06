"""
We want to have flexible arrays, so we do it all in here, and then import from
this module.

This is all probably overly-complicated, and should be removed at first
opportunity to ditch numarray.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Deprecate this.
"""

from config import ytcfg
from logger import ytLogger as ytlog

class ArrayNumTypes:
    def __init__(self):
        pass

u_numpy = True
u_numarray = False

nT = ArrayNumTypes()

myTypes = [ 'Complex128',  'Bool', 'Int32', 'Complex64', 'UInt16', 'Float32',
            'Int64', 'UInt8',  'Int8', 'Complex32',  'UInt64', 'Float64', 
            'UInt32', 'Float128', 'Int16']

import numpy as na
import numpy.linalg as la
import numpy as obj  # Backwards compat
import numpy.core.records as rec
import scipy.ndimage as nd # Moved into scipy
import scipy as sp
import scipy.weave as weave
u_numpy = True
for type in myTypes:
    setattr(nT, type, na.typeNA[type])

# Now define convenience functions

def blankRecordArray(desc, elements):
    """
    Accept a descriptor describing a recordarray, and return one that's full of
    zeros

    This seems like it should be in the numpy distribution...
    """
    blanks = []
    for type in desc['formats']:
        blanks.append(na.zeros(elements, dtype=type))
    return rec.fromarrays(blanks, **desc)

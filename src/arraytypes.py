"""
We want to have flexible arrays, so we do it all in here, and then import from
this module.

This is all probably overly-complicated, and should be removed at first
opportunity to ditch numarray.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Migrate everybody over to numpy, and ditch this
@todo: Ditch the numpy = obj usage
"""

from config import ytcfg
from logger import ytLogger as ytlog

class ArrayNumTypes:
    def __init__(self):
        pass

u_numpy = False
u_numarray = False

nT = ArrayNumTypes()

myTypes = [ 'Complex128',  'Bool', 'Int32', 'Complex64', 'UInt16', 'Float32',
            'Int64', 'UInt8',  'Int8', 'Complex32',  'UInt64', 'Float64', 
            'UInt32', 'Float128', 'Int16']

if not ytcfg.has_option("yt","numarray"):
    ytlog.info("Using NumPy.")
    ytlog.info("Please report problems to mturk@slac.stanford.edu, including a full traceback.")
    ytlog.info("(Check for .flat access to arrays -- it should now be .ravel() !)")
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
else:
    ytlog.warning("Using NumArray.  Many, many things will probably break.  You should use NumPy!")
    import numarray as na
    import numarray.linear_algebra as la
    import numarray.objects as obj
    import numarray.nd_image as nd
    import numarray.records as rec
    u_numarray = True
    for type in myTypes:
        try:
            setattr(nT, type, na.typeDict[type])
        except KeyError:
            pass

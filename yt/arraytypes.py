"""
We want to have flexible arrays, so we do it all in here, and then import from
this module.

This is all probably overly-complicated, and should be removed at first
opportunity to ditch numarray.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import numpy.core.records as rec

# Now define convenience functions

def blankRecordArray(desc, elements):
    """
    Accept a descriptor describing a recordarray, and return one that's full of
    zeros

    This seems like it should be in the numpy distribution...
    """
    blanks = []
    for atype in desc['formats']:
        blanks.append(np.zeros(elements, dtype=atype))
    return rec.fromarrays(blanks, **desc)

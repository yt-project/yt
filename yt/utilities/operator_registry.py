"""
Operation registry class



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import copy
import types

class OperatorRegistry(dict):
    def find(self, op, *args, **kwargs):
        if isinstance(op, types.StringTypes):
            # Lookup, assuming string or hashable object
            op = copy.deepcopy(self[op])
            op.args = args
            op.kwargs = kwargs
        return op

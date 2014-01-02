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

import types

class OperatorRegistry(dict):
    def find(self, op, *args, **kwargs):
        if isinstance(op, types.StringTypes):
            # Lookup, assuming string or hashable object
            op = self[op]
            op.args = args
            op.kwargs = kwargs
        return op

callback_registry = OperatorRegistry()
hf_registry = OperatorRegistry()
quantity_registry = OperatorRegistry()

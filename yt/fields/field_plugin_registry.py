"""
This is a semi-global field plugin registry.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

field_plugins = {}

def register_field_plugin(func):
    name = func.func_name
    if name.startswith("setup_"):
        name = name[len("setup_"):]
    if name.endswith("_fields"):
        name = name[:-len("_fields")]
    field_plugins[name] = func
    # And, we return it, too
    return func

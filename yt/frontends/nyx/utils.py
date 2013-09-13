"""
Utilities for dealing with Nyx data



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

def boxlib_bool_to_int(v):
    try:
        return int(v)
    except ValueError:
        pass
    v = v.upper().strip()
    if v[0] == 'T':
        return 1
    elif v[0] == 'F':
        return 0
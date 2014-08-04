"""
These are particle union objects.  These essentially alias one particle to
another, where the other can be one or several particle types.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import ensure_list

class ParticleUnion(object):
    def __init__(self, name, sub_types):
        self.name = name
        self.sub_types = ensure_list(sub_types)

    def __iter__(self):
        for st in self.sub_types:
            yield st

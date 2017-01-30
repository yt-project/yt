"""
Union structures which can be used to form unions of particles, meshes,
etc. Union is the base class from which trivial named union classes
can be derived



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import ensure_list

class Union(object):
    def __init__(self, name, sub_types):
        self.name = name
        self.sub_types = ensure_list(sub_types)

    def __iter__(self):
        for st in self.sub_types:
            yield st

class MeshUnion(Union):
    def __init__(self, name, sub_types):
        super(MeshUnion, self).__init__(name, sub_types)

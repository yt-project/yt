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

from .unions import Union

class ParticleUnion(Union):
    def __init__(self, name, sub_types):
        super(ParticleUnion, self).__init__(name, sub_types)

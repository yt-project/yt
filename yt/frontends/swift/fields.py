"""
OWLS fields




"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.sph.fields import \
    SPHFieldInfo

class SwiftFieldInfo(SPHFieldInfo):
    _num_neighbors = 48

    def __init__(self, *args, **kwargs):
        super(SwiftFieldInfo,self).__init__( *args, **kwargs )

    def setup_particle_fields(self, ptype):
        """ additional particle fields derived from those in snapshot.
        we also need to add the smoothed fields here b/c setup_fluid_fields
        is called before setup_particle_fields. """

        super(SwiftFieldInfo, self).setup_particle_fields(
            ptype, num_neighbors=self._num_neighbors)

    def setup_fluid_fields(self):
        return

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.gadget.api import GadgetFieldInfo

class ArepoFieldInfo(GadgetFieldInfo):
    known_particle_fields = GadgetFieldInfo.known_particle_fields + \
                            (("smoothing_length", ("code_length", [], None)),)

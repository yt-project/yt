"""
SWIFT fields




"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.sph.fields import \
    SPHFieldInfo

class SwiftFieldInfo(SPHFieldInfo):
    '''
    '''

    def __init__(self, *args, **kwargs):
        super(SwiftFieldInfo,self).__init__( *args, **kwargs )

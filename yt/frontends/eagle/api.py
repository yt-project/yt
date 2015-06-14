"""
API for EAGLE frontend




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .data_structures import \
    EagleDataset, \
    EagleNetworkDataset

from .fields import \
    EagleNetworkFieldInfo

from .io import \
    IOHandlerEagleNetwork

from . import tests

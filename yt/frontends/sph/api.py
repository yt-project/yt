"""
API for SPH frontends




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .eagle.api import \
    EagleDataset, \
    EagleNetworkDataset

from .gadget.api import \
    GadgetDataset

from .http_stream.api import \
    HTTPStreamDataset
    
from .owls.api import \
    OWLSDataset

from .tipsy.api import \
    TipsyDataset

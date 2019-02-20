"""
API for SWIFT frontend

The SWIFT code can be found here:
http://icc.dur.ac.uk/swift/

This frontend is similar in IO to gadget, but has a lot of major differences.

Created by Ashley Kelly (a.j.kelly@durham.ac.uk)

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .data_structures import \
    SwiftDataset

from yt.frontends.sph.fields import \
    SPHFieldInfo

from .io import \
    IOHandlerSwift

from . import tests

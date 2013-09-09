"""
Skeleton-specific IO functions


Authors:
 * Matthew Turk 


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import h5py

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerSkeleton(BaseIOHandler):
    _particle_reader = False
    _data_style = "skeleton"

    def _read_data(self, grid, field):
        # This must return the array, of size/shape grid.ActiveDimensions, that
        # corresponds to 'field'.
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        # If this is not implemented, the IO handler will just slice a
        # _read_data item.
        pass

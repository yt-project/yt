"""
Skeleton-specific IO functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import h5py

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerSkeleton(BaseIOHandler):
    _particle_reader = False
    _data_style = "skeleton"

    def _read_data_set(self, grid, field):
        # This must return the array, of size/shape grid.ActiveDimensions, that
        # corresponds to 'field'.
        pass

    def _read_data_slice(self, grid, field, axis, coord):
        # If this is not implemented, the IO handler will just slice a
        # _read_data_set item.
        pass

"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

import numpy as na

from TransferFunction import TransferFunction, ColorTransferFunction
from yt.amr_utils import PartitionedGrid, VectorPlane, \
                             TransferFunctionProxy
from grid_partitioner import partition_all_grids, partition_grid, \
                             export_partitioned_grids, \
                             import_partitioned_grids
from software_sampler import direct_ray_cast, VolumeRendering
from image_handling import export_rgba, import_rgba, \
                           plot_channel, plot_rgb

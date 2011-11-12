"""
Container file to hold all our Cython routines.  This is to avoid problems with
static linking.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

#cython embedsignature=True
#cython cdivision=True
#cython always_allow_keywords=True

# Set up some imports
import numpy as np
cimport numpy as np
cimport cython

# We include all of our files

include "_amr_utils/DepthFirstOctree.pyx"
include "_amr_utils/Interpolators.pyx"
include "_amr_utils/PointsInVolume.pyx"
include "_amr_utils/RayIntegrators.pyx"
include "_amr_utils/VolumeIntegrator.pyx"
include "_amr_utils/CICDeposit.pyx"
include "_amr_utils/ContourFinding.pyx"
include "_amr_utils/png_writer.pyx"
include "_amr_utils/fortran_reader.pyx"
include "_amr_utils/QuadTree.pyx"
include "_amr_utils/Octree.pyx"
include "_amr_utils/freetype_writer.pyx"
include "_amr_utils/misc_utilities.pyx"
include "_amr_utils/geometry_utils.pyx"

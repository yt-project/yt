"""
Very simple convenience function for importing all the modules, setting up
the namespace and getting the last argument on the command line.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
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


import yt.lagos as lagos
import yt.lagos.hop as hop
import yt.raven as raven
import yt.fido as fido
import numpy as na
import sys

from yt.lagos import EnzoStaticOutput, \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D, \
    add_field, fieldInfo, \
    Clump, write_clump_hierarchy, find_clumps, write_clumps

from yt.raven import PlotCollection, PlotCollectionInteractive, \
    QuiverCallback, ParticleCallback, ContourCallback, \
    GridBoundaryCallback, UnitBoundaryCallback, \
    LinePlotCallback, CuttingQuiverCallback, ClumpContourCallback

try:
    from yt.raven import VolumeRenderingDataCube, \
        VolumeRendering3DProfile, HaloMassesPositionPlot
except ImportError:
    pass

from yt.fido import GrabCollections, OutputCollection

def get_pf():
    return lagos.EnzoStaticOutput(sys.argv[-1])


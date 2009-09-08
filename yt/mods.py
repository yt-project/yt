"""
Very simple convenience function for importing all the modules, setting up
the namespace and getting the last argument on the command line.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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

#
# ALL IMPORTS GO HERE
#

# First module imports
import yt.lagos as lagos
import yt.raven as raven
import yt.fido as fido
import numpy as na
import sys, types
from logger import ytLogger as mylog

# Now individual component imports from lagos
from yt.lagos import EnzoStaticOutput, \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D, \
    add_field, FieldInfo, EnzoFieldInfo, Enzo2DFieldInfo, OrionFieldInfo, \
    Clump, write_clump_hierarchy, find_clumps, write_clumps, \
    OrionStaticOutput, HaloFinder, HOPHaloFinder, FOFHaloFinder

# This is a temporary solution -- in the future, we will allow the user to
# select this via ytcfg.

fieldInfo = EnzoFieldInfo

# Now individual component imports from raven
from yt.raven import PlotCollection, PlotCollectionInteractive, get_multi_plot
from yt.raven.Callbacks import callback_registry
for name, cls in callback_registry.items():
    exec("%s = cls" % name)

# Optional component imports from raven
try:
    from yt.raven import VolumeRenderingDataCube, \
        VolumeRendering3DProfile, HaloMassesPositionPlot
except ImportError:
    pass

import yt.raven.PlotInterface as plots

# Individual imports from Fido
from yt.fido import GrabCollections, OutputCollection

import yt.funcs

from yt.convenience import all_pfs, max_spheres, load

try:
    from yt.sage import *
except ImportError:
    pass

# Some convenience functions to ease our time running scripts
# from the command line

def get_pf():
    return EnzoStaticOutput(sys.argv[-1])

def get_pc():
    return PlotCollection(EnzoStaticOutput(sys.argv[-1]))

mylog.info("Welcome to YT.")

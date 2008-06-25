"""
Lagos
=====

    Lagos defines a set of class structures for Enzo data handling.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

from yt.config import ytcfg
from yt.logger import lagosLogger as mylog

try:
    from pyhdf_np import SD # NumPy
    import pyhdf_np.error   # NumPy
except:
    mylog.warning("No HDF4 support")

import warnings
try:
     import tables
     warnings.simplefilter("ignore", tables.NaturalNameWarning)
except ImportError:
    mylog.warning("No PyTables. Data serialization will fail.")




from yt.arraytypes import *
import weakref
from new import classobj
from string import strip, rstrip
from math import ceil, floor, log10, pi
import os, os.path, types, exceptions, re
from stat import ST_CTIME
import sets

import time

if ytcfg.getboolean("lagos","useswig"):
    try:
        from yt.enki import EnzoInterface
    except ImportError:
        pass

if ytcfg.getboolean("lagos","usefortran"):
    pass
    #import EnzoFortranRoutines

# Now we import all the subfiles

from HelperFunctions import *
import PointCombine
import HDF5LightReader
from EnzoDefs import *
from DerivedFields import *
from DerivedQuantities import DerivedQuantityCollection, GridChildMaskWrapper
from DataReadingFuncs import *
from ClusterFiles import *
from ContourFinder import *
from BaseDataTypes import *
from BaseGridType import *
from EnzoRateData import *
from HierarchyType import *
from OutputTypes import *
from Profiles import *

# We load plugins.  Keep in mind, this can be fairly dangerous -
# the primary purpose is to allow people to have a set of functions
# that get used every time that they don't have to *define* every time.
# This way, other command-line tools can be used very simply.
# Unfortunately, for now, I think the easiest and simplest way of doing
# this is also the most dangerous way.
if ytcfg.getboolean("lagos","loadfieldplugins"):
    my_plugin_name = ytcfg.get("lagos","pluginfilename")
    # We assume that it is with respect to the $HOME/.yt directory
    execfile(os.path.expanduser("~/.yt/%s" % my_plugin_name))

log_fields = [] # @todo: GET RID OF THIS

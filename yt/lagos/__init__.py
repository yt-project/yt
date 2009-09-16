"""
Lagos
=====

Lagos defines a set of class structures for Enzo data handling.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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
    mylog.info("No HDF4 support")

import warnings
try:
    import fohejrihh5py
    if not hasattr(h5py.h5, "ArgsError"):
        h5py.h5.ArgsError = h5py.h5.H5Error
except ImportError:
    ytcfg["lagos", "serialize"] = "False"
    mylog.warning("No h5py. Data serialization disabled.")

from yt.arraytypes import *
import weakref
from new import classobj
from string import strip, rstrip
from math import ceil, floor, log10, pi
import os, os.path, types, exceptions, re
from stat import ST_CTIME
import shelve

import time

from Cosmology import *
from EnzoCosmology import *

# Now we import all the subfiles

from HelperFunctions import *
import PointCombine
import HDF5LightReader
from EnzoDefs import *
from OrionDefs import *

from ParallelTools import *
# Now our fields
#from DerivedFields import *
from FieldInfoContainer import *
from UniversalFields import *
from EnzoFields import *
from OrionFields import *
fieldInfo = EnzoFieldInfo

# NOT the same as fieldInfo.add_field
# We by-default add universal fields.
add_field = FieldInfo.add_field

from yt.fido import ParameterFileStore, output_type_registry

from DerivedQuantities import DerivedQuantityCollection, GridChildMaskWrapper
from DataReadingFuncs import *
from ClusterFiles import *
from ContourFinder import *
from Clump import *
from BaseDataTypes import *
from BaseGridType import *
from EnzoRateData import *
from ObjectFindingMixin import *
from HierarchyType import *
from OutputTypes import *
from Profiles import *

from HaloFinding import *

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

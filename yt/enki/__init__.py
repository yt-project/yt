"""
Enki
====

Enki is the package used to create data, and instantiate runs.  It supports
creating Enzo Problems, and then using SWIG-interfaced Enzo calls to actually
create the data for those problems.  Additionally, facilities are being
developed to use Enki to directly execute runs.

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

import sys
from yt.logger import enkiLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *

# Now we import the SWIG enzo interface
# Note that we're going to try super-hard to get the one that's local for the
# user

sp = sys.path

if ytcfg.has_option("SWIG", "EnzoInterfacePath"):
    swig_path = ytcfg.get("SWIG","EnzoInterfacePath")
    mylog.info("Using %s as path to SWIG Interface", swig_path)
    sys.path = sys.path[:1] + [swig_path] + sys.path[1:] # We want '' to be the first

try:
    import EnzoInterface
    mylog.debug("Imported EnzoInterface successfully")
    has_SWIG = True
except ImportError, e:
    mylog.warning("EnzoInterface failed to import; all SWIG actions will fail")
    mylog.warning("(%s)", e)
    has_SWIG = False

sys.path = sp

from EnkiDefs import *

from HostProfiles import *
from CreateNewProblem import *
import mes
from Recompile import *

#from yt.lagos import * # Do we need this?

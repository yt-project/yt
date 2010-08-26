"""
API for yt.visualization

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: J.S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

import time, types, string, os

from yt.funcs import *

from yt.config import ytcfg
from yt.utilities.logger import ravenLogger as mylog
#from yt.arraytypes import *

vm_axis_names = {0:'x', 1:'y', 2:'z', 3:'dx', 4:'dy'}

import plot_types
be = plot_types

#from plot_modifications import *
#from fixed_resolution import *

#from plot_collection import *

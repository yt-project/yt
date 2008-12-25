"""
Some convenience functions, objects, and iterators

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

import glob

# Named imports
import yt.lagos as lagos
import yt.raven as raven
from yt.funcs import *
import numpy as na
import os.path, inspect, types
from functools import wraps
from yt.logger import ytLogger as mylog

def all_pfs(max_depth=1, name_spec="*.hierarchy"):
    list_of_names = []
    for i in range(max_depth):
        bb = list('*' * i) + [name_spec]
        list_of_names += glob.glob(os.path.join(*bb))
    list_of_names.sort(key=lambda b: os.path.basename(b))
    for fn in list_of_names:
        yield lagos.EnzoStaticOutput(fn[:-10])

def max_spheres(width, unit, **kwargs):
    for pf in all_pfs(**kwargs):
        v, c = pf.h.find_max("Density")
        yield pf.h.sphere(c, width/pf[unit])

"""
Some convenience functions, objects, and iterators

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

import glob

# Named imports
import yt.lagos as lagos
import yt.raven as raven
from yt.funcs import *
import numpy as na
import os.path, inspect, types
from functools import wraps
from yt.logger import ytLogger as mylog
from yt.fido import output_type_registry

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

def load(*args ,**kwargs):
    candidates = []
    for n, c in output_type_registry.items():
        if n is None: continue
        if c._is_valid(*args, **kwargs): candidates.append(n)
    if len(candidates) == 1:
        return output_type_registry[candidates[0]](*args, **kwargs)
    if len(candidates) == 0:
        mylog.error("Couldn't figure out output type for %s", args[0])
        return None
    mylog.error("Multiple output type candidates for %s:", args[0])
    for c in candidates:
        mylog.error("    Possible: %s", c)
    return None

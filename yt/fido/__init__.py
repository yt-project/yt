"""
Fido is a module for moving and storing data.  It is designed to watch for
outputs, as well as transparent archiving them during and after a simulation
run.  Additionally, function handlers can be called to deal with output upon
its creation.

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

from yt.logger import fidoLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *

import os, os.path, shutil, time, sys, glob, types

OTHERFILES=["rates.out","cool_rates.out"]
NEW_OUTPUT_CREATED = ytcfg.get("fido","NewOutputCreated")
GlobPatterns = ytcfg.get("fido","GlobPatterns").split(",")
NewDirectoryPattern = ytcfg.get("fido","NewDirectoryPattern",raw=True)

# Let's define some useful functions.  Maybe these should go elsewhere?

def get_parameter_line(filename, parameter):
    f = open(filename)
    lines = filter(lambda a: a.startswith(parameter), f)
    if len(lines) > 1:
        raise KeyError, "More than one line matches that parameter!"
    else: return lines[0]

def get_parent_dir(filename):
    return os.path.normpath( \
            os.path.split(os.path.dirname(os.path.abspath(filename)))[0])

from OutputCollection import *
from FileHandling import *
from OutputWatcher import *
from RunStandalones import *

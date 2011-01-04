"""
Various definitions for various other modules and routines - modified
from Orion frontend.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Chris Malone <chris.m.malone@gmail.com>
Affiliation: SUNY Stony Brook
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-20010 J.S. Oishi.  All Rights Reserved.

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
from yt.funcs import *

parameterTypes = {"dimensionality": int,
                 "refine_by": int
                 }

# converts the Maestro inputs file name to the Enzo/yt name expected
# throughout the code. key is Maestro name, value is Enzo/yt equivalent
maestro2enzoDict = {"dm_in": "dimensionality",
                    "ref_ratio": "refine_by"
                  }

yt2maestroFieldsDict = {}
maestro2ytFieldsDict = {}

maestro_FAB_header_pattern = r"^FAB \(\((\d+), \([0-9 ]+\)\),\(\d+, \(([0-9 ]+)\)\)\)\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n"

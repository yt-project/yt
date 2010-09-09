"""
A mechanism for listing available analysis modules.

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

import os

def get_available_modules():
    modpath = os.path.abspath(os.path.dirname(__file__))
    available_modules = []
    for d in [os.path.join(modpath, f) for f in os.listdir(modpath)]:
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, "api.py")):
            available_modules.append(os.path.basename(d))
    return available_modules

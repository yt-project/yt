"""
Utilities for dealing with Nyx data

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Casey W. Stark, Matthew Turk.  All Rights Reserved.

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

def boxlib_bool_to_int(v):
    try:
        return int(v)
    except ValueError:
        pass
    v = v.upper().strip()
    if v[0] == 'T':
        return 1
    elif v[0] == 'F':
        return 0
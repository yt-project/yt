"""
Chombo-specific fields

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: UC Berkeley
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 J. S. Oishi, Matthew Turk.  All Rights Reserved.

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

from UniversalFields import *
class ChomboFieldContainer(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
ChomboFieldInfo = ChomboFieldContainer()
add_chombo_field = ChomboFieldInfo.add_field

add_field = add_chombo_field

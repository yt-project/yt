"""
RAMSES-specific fields

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
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

from yt.data_objects.field_info_container import \
    CodeFieldInfoContainer, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields

class RAMSESFieldContainer(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
RAMSESFieldInfo = RAMSESFieldContainer()
add_ramses_field = RAMSESFieldInfo.add_field

add_field = add_ramses_field

known_ramses_fields = [
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
]

for f in known_ramses_fields:
    if f not in RAMSESFieldInfo:
        add_field(f, function=lambda a,b: None, take_log=True,
                  validators = [ValidateDataField(f)])

def _convertDensity(data):
    return data.convert("Density")
RAMSESFieldInfo["Density"]._units = r"\rm{g}/\rm{cm}^3"
RAMSESFieldInfo["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
RAMSESFieldInfo["Density"]._convert_function=_convertDensity

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = RAMSESFieldInfo["%s-velocity" % ax]
    f._units = r"\rm{cm}/\rm{s}"
    f._convert_function = _convertVelocity
    f.take_log = False


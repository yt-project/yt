"""
Maestro-specific fields - borrows heavily from Orion frontend.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: UC Berkeley
Author: Chris Malone <chris.m.malone@gmail.com>
Affiliation: SUNY Stony Brook
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 J. S. Oishi, Matthew Turk.  All Rights Reserved.

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
from yt.utilities.physical_constants import \
    mh, kboltz
from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields

KnownMaestroFields = FieldInfoContainer()
add_maestro_field = KnownMaestroFields.add_field

MaestroFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = MaestroFieldInfo.add_field

add_field("density", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("density")],
          units="g/cm**3")

translation_dict = {
    "x-velocity": "x_vel",
    "y-velocity": "y_vel",
    "z-velocity": "z_vel",
    "Density":    "density",
    "Temperature": "tfromp"
}

def _generate_translation(mine, theirs):
    add_field(theirs, function=lambda a, b: b[mine], take_log=True)

for f, v in translation_dict.items():
    if v not in MaestroFieldInfo:
        add_field(v, function=lambda a,b: None, take_log=False,
                  validators = [ValidateDataField(v)])
#    print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)


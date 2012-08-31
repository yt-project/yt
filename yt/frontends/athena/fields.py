"""
Athena-specific fields

Author: Samuel W. Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Samuel W. Skillman, Matthew Turk, J. S. Oishi.  
  All Rights Reserved.

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
    FieldInfoContainer, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType, \
    NullFunc, \
    TranslationFunc
import yt.data_objects.universal_fields

log_translation_dict = {"Density": "density",
                        "Pressure": "pressure"}

translation_dict = {"x-velocity": "velocity_x",
                    "y-velocity": "velocity_y",
                    "z-velocity": "velocity_z"}
                    
# translation_dict = {"mag_field_x": "cell_centered_B_x ",
#                     "mag_field_y": "cell_centered_B_y ",
#                     "mag_field_z": "cell_centered_B_z "}

AthenaFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = AthenaFieldInfo.add_field

KnownAthenaFields = FieldInfoContainer()
add_athena_field = KnownAthenaFields.add_field

add_athena_field("density", function=NullFunc, take_log=True,
          units=r"\rm{g}/\rm{cm}^3",
          projected_units =r"\rm{g}/\rm{cm}^2")

add_athena_field("specific_energy", function=NullFunc, take_log=True,
          units=r"\rm{erg}/\rm{g}")

add_athena_field("pressure", function=NullFunc, take_log=True,
          units=r"\rm{erg}/\rm{g}")

add_athena_field("velocity_x", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_athena_field("velocity_y", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_athena_field("velocity_z", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_athena_field("mag_field_x", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_athena_field("mag_field_y", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_athena_field("mag_field_z", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

for f,v in log_translation_dict.items():
    add_field(f, TranslationFunc(v), take_log=True)

for f,v in translation_dict.items():
    add_field(f, TranslationFunc(v), take_log=False)


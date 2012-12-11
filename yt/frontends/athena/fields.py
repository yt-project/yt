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

import numpy as np
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
from yt.utilities.physical_constants import \
    kboltz,mh
import yt.data_objects.universal_fields

log_translation_dict = {}

translation_dict = {"Density": "density"}

AthenaFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = AthenaFieldInfo.add_field

KnownAthenaFields = FieldInfoContainer()
add_athena_field = KnownAthenaFields.add_field

add_athena_field("density", function=NullFunc, take_log=False,
                 units=r"", projected_units =r"")

add_athena_field("pressure", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("velocity_x", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("velocity_y", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("velocity_z", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("momentum_x", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("momentum_y", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("momentum_z", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("cell_centered_B_x", function=NullFunc, take_log=False,
                 units=r"", display_name=r"$\rm{cell\ centered\ B_x}$")

add_athena_field("cell_centered_B_y", function=NullFunc, take_log=False,
                 units=r"", display_name=r"$\rm{cell\ centered\ B_y}$")

add_athena_field("cell_centered_B_z", function=NullFunc, take_log=False,
                 units=r"", display_name=r"$\rm{cell\ centered\ B_z}$")

for f,v in log_translation_dict.items():
    add_field(f, TranslationFunc(v), take_log=True)

for f,v in translation_dict.items():
    add_field(f, TranslationFunc(v), take_log=False)

# In Athena, conservative or primitive variables may be written out.
# By default, yt concerns itself with primitive variables. The following
# field definitions allow for conversions to primitive variables in the
# case that the file contains the conservative ones.

def _xvelocity(field, data) :
    try:
        xvel = data["velocity_x"]
    except:
        xvel = data["momentum_x"]/data["Density"]
    return xvel
add_field("x-velocity", function=_xvelocity, take_log=False,
          units=r"")

def _yvelocity(field, data) :
    try:
        yvel = data["velocity_y"]
    except:
        yvel = data["momentum_y"]/data["Density"]
    return yvel
add_field("y-velocity", function=_yvelocity, take_log=False,
          units=r"")

def _zvelocity(field, data) :
    try:
        zvel = data["velocity_z"]
    except:
        zvel = data["momentum_z"]/data["Density"]
    return zvel
add_field("z-velocity", function=_zvelocity, take_log=False,
          units=r"")

def _Temperature(field, data):
    if data.has_field_parameter("mu") :
        mu = data.get_field_parameter("mu")
    else:
        mu = 0.6
    return mu*mh*data["Pressure"]/data["Density"]/kboltz
add_field("Temperature", function=_Temperature, take_log=False,
          units=r"\rm{K}")

def _Bx(field, data):
    factor = np.sqrt(4.*np.pi)
    return data['cell_centered_B_x']*factor
add_field("Bx", function=_Bx, take_log=False,
          units=r"\rm{Gauss}", display_name=r"B_x")

def _By(field, data):
    factor = np.sqrt(4.*np.pi)
    return data['cell_centered_B_y']*factor
add_field("By", function=_By, take_log=False,
          units=r"\rm{Gauss}", display_name=r"B_y")

def _Bz(field, data):
    factor = np.sqrt(4.*np.pi)
    return data['cell_centered_B_z']*factor
add_field("Bz", function=_Bz, take_log=False,
          units=r"\rm{Gauss}", display_name=r"B_z")


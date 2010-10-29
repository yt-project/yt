"""
ART-specific fields

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
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, mass_hydrogen_cgs

class ARTFieldContainer(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
ARTFieldInfo = ARTFieldContainer()
add_art_field = ARTFieldInfo.add_field

add_field = add_art_field

translation_dict = {"Density":"density",
                    "x-velocity":"velocity_x",
                    "y-velocity":"velocity_y",
                    "z-velocity":"velocity_z",
                    "Pressure":"pressure",
                    "Metallicity":"metallicity",
                    "Gas_Energy":"Gas_Energy",
                    "Total_Energy":"Total_Energy",
                   }

def _generate_translation(mine, theirs):
    add_field(theirs, function=lambda a, b: b[mine], take_log=True)

for f,v in translation_dict.items():
    if v not in ARTFieldInfo:
        add_field(v, function=lambda a,b: None, take_log=False,
                  validators = [ValidateDataField(v)])
    #print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)

def _convertDensity(data):
    return data.convert("Density")
ARTFieldInfo["Density"]._units = r"\rm{g}/\rm{cm}^3"
ARTFieldInfo["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
ARTFieldInfo["Density"]._convert_function=_convertDensity

def _convertEnergy(data):
    return data.convert("Gas_Energy")
ARTFieldInfo["Gas_Energy"]._units = r"\rm{ergs}/\rm{g}"
ARTFieldInfo["Gas_Energy"]._convert_function=_convertEnergy

def _Temperature(field, data):
    tr = (data["Gamma"] - 1.0) * data["Gas_Energy"]
    tr *= mass_hydrogen_cgs * data.pf["wmu"] / boltzmann_constant_cgs
    return tr
add_field("Temperature", function=_Temperature, units = r"\mathrm{K}")

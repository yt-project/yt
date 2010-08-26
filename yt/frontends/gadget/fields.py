"""
Gadget-specific fields

Author: Christopher E Moody <juxtaposicion@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Christopher E Moody, Matthew Turk.  All Rights Reserved.

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

class GadgetFieldContainer(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
GadgetFieldInfo = GadgetFieldContainer()
add_gadget_field = GadgetFieldInfo.add_field

add_field = add_gadget_field

translation_dict = {"Position": "POS",
                    "Velocity": "VEL",
                    "ID": "ID",
                    "Mass":"MASS"
                   }

#for f,v in translation_dict.items():
#    add_field(f, function=lambda a,b: None, take_log=True,
#        validators = [ValidateDataField(v)],
#        units=r"\rm{cm}")
#    add_field(v, function=lambda a,b: None, take_log=True,
#        validators = [ValidateDataField(v)],
#        units=r"\rm{cm}")
          
          
add_field("Density", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("POS")],
          units=r"\rm{cm}")

add_field("VEL", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("VEL")],
          units=r"")

add_field("ID", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("ID")],
          units=r"")

add_field("MASS", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("MASS")],
          units=r"\rm{g}")

add_field("U", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("U")],
          units=r"")

add_field("NE", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("NE")],
          units=r"")

add_field("POT", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("POT")],
          units=r"")

add_field("ACCE", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("ACCE")],
          units=r"")

add_field("ENDT", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("ENDT")],
          units=r"")

add_field("TSTP", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("TSTP")],
          units=r"")


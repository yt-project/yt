"""MOAB-specific fields


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

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

MoabFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = MoabFieldInfo.add_field

KnownMoabFields = FieldInfoContainer()
add_moab_field = KnownMoabFields.add_field

add_moab_field("density", function=NullFunc, take_log=True,
          units=r"\rm{g}/\rm{cm}^3",
          projected_units =r"\rm{g}/\rm{cm}^2")

add_moab_field("specific_energy", function=NullFunc, take_log=True,
          units=r"\rm{erg}/\rm{g}")

add_moab_field("pressure", function=NullFunc, take_log=True,
          units=r"\rm{erg}/\rm{g}")

add_moab_field("velocity_x", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_moab_field("velocity_y", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_moab_field("velocity_z", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_moab_field("mag_field_x", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_moab_field("mag_field_y", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

add_moab_field("mag_field_z", function=NullFunc, take_log=False,
          units=r"\rm{cm}/\rm{s}")

for f,v in log_translation_dict.items():
    add_field(f, TranslationFunc(v), take_log=True)

for f,v in translation_dict.items():
    add_field(f, TranslationFunc(v), take_log=False)


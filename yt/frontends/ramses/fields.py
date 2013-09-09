"""
RAMSES-specific fields


Authors:
 * Matthew Turk 


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
    ValidateGridType
import yt.data_objects.universal_fields


KnownRAMSESFields = FieldInfoContainer()
add_ramses_field = KnownRAMSESFields.add_field

RAMSESFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = RAMSESFieldInfo.add_field

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


"""
Maestro-specific fields - borrows heavily from Orion frontend.


Authors:
 * J. S. Oishi 
 * Chris Malone 


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
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
          units=r"\rm{g}/\rm{cm}^3")
MaestroFieldInfo["density"]._projected_units =r"\rm{g}/\rm{cm}^2"

translation_dict = {"x-velocity": "x_vel",
                    "y-velocity": "y_vel",
                    "z-velocity": "z_vel",
                    "Density": "density",
                    "Temperature": "tfromp"
                   }

def _generate_translation(mine, theirs):
    add_field(theirs, function=lambda a, b: b[mine], take_log=True)

for f,v in translation_dict.items():
    if v not in MaestroFieldInfo:
        add_field(v, function=lambda a,b: None, take_log=False,
                  validators = [ValidateDataField(v)])
#    print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)


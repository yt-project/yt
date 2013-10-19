"""
Tiger-specific fields



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

KnownTigerFields = FieldInfoContainer()
add_tiger_field = KnownTigerFields.add_field

TigerFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = TigerFieldInfo.add_field


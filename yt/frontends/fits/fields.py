"""
FITS-specific fields
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.utilities.exceptions import *
from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    NullFunc, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.fields.universal_fields
KnownFITSFields = FieldInfoContainer()
add_fits_field = KnownFITSFields.add_field

FITSFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = FITSFieldInfo.add_field


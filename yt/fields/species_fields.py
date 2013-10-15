"""
Fields based on species of molecules or atoms.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

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
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions, \
    standard_particle_fields
from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs
from yt.funcs import *

import yt.utilities.lib as amr_utils

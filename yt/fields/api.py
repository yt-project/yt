""" 
API for yt.fields


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .field_plugin_registry import \
    register_field_plugin, \
    field_plugins

from . import angular_momentum
from . import astro_fields
from . import cosmology_fields
from . import fluid_fields
from . import fluid_vector_fields
from . import magnetic_field
from . import geometric_fields
from . import particle_fields
#from . import species_fields
from . import vector_operations
from . import local_fields
from . import my_plugin_fields

from .local_fields import add_field, derived_field


from .derived_field import \
    DerivedField, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
from .field_detector import \
    FieldDetector
from .field_info_container import \
    FieldInfoContainer

"""
Field-related exceptions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


class ValidationException(Exception):
    pass

class NeedsGridType(ValidationException):
    def __init__(self, ghost_zones = 0, fields=None):
        self.ghost_zones = ghost_zones
        self.fields = fields
    def __str__(self):
        return "(%s, %s)" % (self.ghost_zones, self.fields)

class NeedsOriginalGrid(NeedsGridType):
    def __init__(self):
        self.ghost_zones = 0

class NeedsDataField(ValidationException):
    def __init__(self, missing_fields):
        self.missing_fields = missing_fields
    def __str__(self):
        return "(%s)" % (self.missing_fields)

class NeedsProperty(ValidationException):
    def __init__(self, missing_properties):
        self.missing_properties = missing_properties
    def __str__(self):
        return "(%s)" % (self.missing_properties)

class NeedsParameter(ValidationException):
    def __init__(self, missing_parameters):
        self.missing_parameters = missing_parameters
    def __str__(self):
        return "(%s)" % (self.missing_parameters)

class NeedsParameterValue(ValidationException):
    def __init__(self, parameter_values):
        self.parameter_values = parameter_values

class NeedsConfiguration(ValidationException):
    def __init__(self, parameter, value):
        self.parameter = parameter
        self.value = value
    def __str__(self):
        return "(Needs %s = %s)" % (self.parameter, self.value)

class FieldUnitsError(Exception):
    pass


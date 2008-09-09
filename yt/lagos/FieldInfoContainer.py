"""
The basic field info container resides here.  These classes, code specific and
universal, are the means by which we access fields across YT, both derived and
native.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

import types
import numpy as na
import inspect
import copy
import itertools

from yt.funcs import *

class FieldInfoContainer(object): # We are all Borg.
    _shared_state = {}
    _universal_field_list = {}
    def __new__(cls, *args, **kwargs):
        self = object.__new__(cls, *args, **kwargs)
        self.__dict__ = cls._shared_state
        return self
    def __getitem__(self, key):
        if key in self._universal_field_list:
            return self._universal_field_list[key]
        raise KeyError
    def keys(self):
        return self._universal_field_list.keys()
    def __iter__(self):
        return self._universal_field_list.iterkeys()
    def __setitem__(self, key, val):
        self._universal_field_list[key] = val
    def add_field(self, name, function = None, **kwargs):
        if function == None:
            if kwargs.has_key("function"):
                function = kwargs.pop("function")
            else:
                # This will fail if it does not exist,
                # which is our desired behavior
                function = eval("_%s" % name)
        self[name] = DerivedField(name, function, **kwargs)
FieldInfo = FieldInfoContainer()
add_field = FieldInfo.add_field

class CodeFieldInfoContainer(FieldInfoContainer):
    def __setitem__(self, key, val):
        self._field_list[key] = val
    def __iter__(self):
        return itertools.chain(self._field_list.iterkeys(),
                        self._universal_field_list.iterkeys())
    def keys(self):
        return set(self._field_list.keys() + self._universal_field_list.keys())
    def __getitem__(self, key):
        if key in self._field_list:
            return self._field_list[key]
        if key in self._universal_field_list:
            return self._universal_field_list[key]
        raise KeyError

class EnzoFields(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
EnzoFieldInfo = EnzoFields()
add_enzo_field = EnzoFieldInfo.add_field

class OrionFields(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
OrionFieldInfo = OrionFields()
add_orion_field = OrionFieldInfo.add_field

class ValidationException(Exception):
    pass

class NeedsGridType(ValidationException):
    def __init__(self, ghost_zones = 0, fields=None):
        self.ghost_zones = ghost_zones
        self.fields = fields

class NeedsDataField(ValidationException):
    def __init__(self, missing_fields):
        self.missing_fields = missing_fields

class NeedsProperty(ValidationException):
    def __init__(self, missing_properties):
        self.missing_properties = missing_properties

class NeedsParameter(ValidationException):
    def __init__(self, missing_parameters):
        self.missing_parameters = missing_parameters

class FieldDetector(defaultdict):
    Level = 1
    NumberOfParticles = 0
    def __init__(self, nd = 16, pf = None):
        self.nd = nd
        self.ActiveDimensions = [nd,nd,nd]
        self.LeftEdge = [0.0,0.0,0.0]
        self.RightEdge = [1.0,1.0,1.0]
        self.dx = self.dy = self.dz = na.array([1.0])
        self.fields = []
        if pf is None:
            pf = defaultdict(lambda: 1)
        self.pf = pf
        self.requested = []
        self.requested_parameters = []
        defaultdict.__init__(self, lambda: na.ones((nd,nd,nd)))
    def __missing__(self, item):
        if FieldInfo.has_key(item) and \
            FieldInfo[item]._function.func_name != '<lambda>':
            try:
                vv = FieldInfo[item](self)
            except NeedsGridType, exc:
                ngz = exc.ghost_zones
                nfd = FieldDetector(self.nd+ngz*2)
                vv = FieldInfo[item](nfd)[ngz:-ngz,ngz:-ngz,ngz:-ngz]
                for i in vv.requested:
                    if i not in self.requested: self.requested.append(i)
                for i in vv.requested_parameters:
                    if i not in self.requested_parameters: self.requested_parameters.append(i)
            if vv is not None:
                self[item] = vv
                return self[item]
        self.requested.append(item)
        return defaultdict.__missing__(self, item)
    def get_field_parameter(self, param):
        self.requested_parameters.append(param)
        if param in ['bulk_velocity','center','height_vector']:
            return na.array([0,0,0])
        else:
            return 0.0
    _spatial = True
    _num_ghost_zones = 0
    id = 1
    def has_field_parameter(self, param): return True
    def convert(self, item): return 1

class DerivedField:
    def __init__(self, name, function,
                 convert_function = None,
                 units = "", projected_units = "",
                 take_log = True, validators = None,
                 particle_type = False, vector_field=False,
                 display_field = True, not_in_all=False,
                 projection_conversion = "cm"):
        """
        This is the base class used to describe a cell-by-cell derived field.

        :param name: is the name of the field.
        :param function: is a function handle that defines the field
        :param convert_function: must convert to CGS, if it needs to be done
        :param units: is a mathtext-formatted string that describes the field
        :param projected_units: if we display a projection, what should the units be?
        :param take_log: describes whether the field should be logged
        :param validators: is a list of :class:`FieldValidator` objects
        :param particle_type: is this field based on particles?
        :param vector_field: describes the dimensionality of the field
        :param display_field: governs its appearance in the dropdowns in reason
        :param not_in_all: is used for baryon fields from the data that are not in
                           all the grids
        :param projection_conversion: which unit should we multiply by in a
                                      projection?
        """
        self.name = name
        self._function = function
        if validators:
            self.validators = ensure_list(validators)
        else:
            self.validators = []
        self.take_log = take_log
        self._units = units
        self._projected_units = projected_units
        if not convert_function:
            convert_function = lambda a: 1.0
        self._convert_function = convert_function
        self.particle_type = particle_type
        self.vector_field = vector_field
        self.projection_conversion = projection_conversion
        self.display_field = display_field
        self.not_in_all = not_in_all
    def check_available(self, data):
        for validator in self.validators:
            validator(data)
        # If we don't get an exception, we're good to go
        return True
    def get_dependencies(self, *args, **kwargs):
        e = FieldDetector(*args, **kwargs)
        if self._function.func_name == '<lambda>':
            e.requested.append(self.name)
        else:
            self(e)
        return e
    def get_units(self):
        return self._units
    def get_projected_units(self):
        return self._projected_units
    def __call__(self, data):
        ii = self.check_available(data)
        original_fields = data.fields[:] # Copy
        dd = self._function(self, data)
        dd *= self._convert_function(data)
        for field_name in data.fields:
            if field_name not in original_fields:
                del data[field_name]
        return dd
    def get_source(self):
        return inspect.getsource(self._function)

class FieldValidator(object):
    pass

class ValidateParameter(FieldValidator):
    def __init__(self, parameters):
        FieldValidator.__init__(self)
        self.parameters = ensure_list(parameters)
    def __call__(self, data):
        doesnt_have = []
        for p in self.parameters:
            if not data.has_field_parameter(p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsParameter(doesnt_have)
        return True

class ValidateDataField(FieldValidator):
    def __init__(self, field):
        FieldValidator.__init__(self)
        self.fields = ensure_list(field)
    def __call__(self, data):
        doesnt_have = []
        if isinstance(data, FieldDetector): return True
        for f in self.fields:
            if f not in data.hierarchy.field_list:
                doesnt_have.append(f)
        if len(doesnt_have) > 0:
            raise NeedsDataField(doesnt_have)
        return True

class ValidateProperty(FieldValidator):
    def __init__(self, prop):
        FieldValidator.__init__(self)
        self.prop = ensure_list(prop)
    def __call__(self, data):
        doesnt_have = []
        for p in self.prop:
            if not hasattr(data,p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsProperty(doesnt_have)
        return True

class ValidateSpatial(FieldValidator):
    def __init__(self, ghost_zones = 0, fields=None):
        FieldValidator.__init__(self)
        self.ghost_zones = ghost_zones
        self.fields = fields
    def __call__(self, data):
        # When we say spatial information, we really mean
        # that it has a three-dimensional data structure
        if isinstance(data, FieldDetector): return True
        if not data._spatial:
            raise NeedsGridType(self.ghost_zones,self.fields)
        if self.ghost_zones == data._num_ghost_zones:
            return True
        raise NeedsGridType(self.ghost_zones)

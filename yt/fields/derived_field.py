"""
Derived field base class.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import contextlib
import inspect

from yt.extern.six import string_types, PY2
from yt.funcs import \
    ensure_list
from .field_exceptions import \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter, \
    NeedsParameterValue, \
    FieldUnitsError
from .field_detector import \
    FieldDetector
from yt.units.unit_object import \
    Unit
import yt.units.dimensions as ytdims
from yt.utilities.exceptions import \
    YTFieldNotFound


def TranslationFunc(field_name):
    def _TranslationFunc(field, data):
        # We do a bunch of in-place modifications, so we will copy this.
        return data[field_name].copy()
    _TranslationFunc.alias_name = field_name
    return _TranslationFunc

def NullFunc(field, data):
    raise YTFieldNotFound(field.name)

class DerivedField(object):
    """
    This is the base class used to describe a cell-by-cell derived field.

    Parameters
    ----------

    name : str
       is the name of the field.
    function : callable
       A function handle that defines the field.  Should accept
       arguments (field, data)
    units : str
       A plain text string encoding the unit, or a query to a unit system of
       a dataset. Powers must be in python syntax (** instead of ^). If set
       to "auto" the units will be inferred from the units of the return
       value of the field function, and the dimensions keyword must also be
       set (see below).
    take_log : bool
       Describes whether the field should be logged
    validators : list
       A list of :class:`FieldValidator` objects
    particle_type : bool
       Is this a particle (1D) field?
    vector_field : bool
       Describes the dimensionality of the field.  Currently unused.
    display_field : bool
       Governs its appearance in the dropdowns in Reason
    not_in_all : bool
       Used for baryon fields from the data that are not in all the grids
    display_name : str
       A name used in the plots
    output_units : str
       For fields that exist on disk, which we may want to convert to other
       fields or that get aliased to themselves, we can specify a different
       desired output unit than the unit found on disk.
    dimensions : str or object from yt.units.dimensions
       The dimensions of the field, only needed if units="auto" and only used
       for error checking.
    """
    def __init__(self, name, function, units=None,
                 take_log=True, validators=None,
                 particle_type=False, vector_field=False, display_field=True,
                 not_in_all=False, display_name=None, output_units=None,
                 dimensions=None, ds=None):
        self.name = name
        self.take_log = take_log
        self.display_name = display_name
        self.not_in_all = not_in_all
        self.display_field = display_field
        self.particle_type = particle_type
        self.vector_field = vector_field
        self.ds = ds

        self._function = function

        if validators:
            self.validators = ensure_list(validators)
        else:
            self.validators = []

        # handle units
        if units is None:
            self.units = ''
        elif isinstance(units, string_types):
            if units.lower() == 'auto':
                if dimensions is None:
                    raise RuntimeError("To set units='auto', please specify the dimensions "
                                       "of the field with dimensions=<dimensions of field>!")
                self.units = None
            else:
                self.units = units
        elif isinstance(units, Unit):
            self.units = str(units)
        else:
            raise FieldUnitsError("Cannot handle units '%s' (type %s)." \
                                  "Please provide a string or Unit " \
                                  "object." % (units, type(units)) )
        if output_units is None:
            output_units = self.units
        self.output_units = output_units

        if isinstance(dimensions, string_types):
            dimensions = getattr(ytdims, dimensions)
        self.dimensions = dimensions

    def _copy_def(self):
        dd = {}
        dd['name'] = self.name
        dd['units'] = self.units
        dd['take_log'] = self.take_log
        dd['validators'] = list(self.validators)
        dd['particle_type'] = self.particle_type
        dd['vector_field'] = self.vector_field
        dd['display_field'] = True
        dd['not_in_all'] = self.not_in_all
        dd['display_name'] = self.display_name
        return dd

    def get_units(self):
        if self.ds is not None:
            u = Unit(self.units, registry=self.ds.unit_registry)
        else:
            u = Unit(self.units)
        return u.latex_representation()

    def get_projected_units(self):
        if self.ds is not None:
            u = Unit(self.units, registry=self.ds.unit_registry)
        else:
            u = Unit(self.units)
        return (u*Unit('cm')).latex_representation()

    def check_available(self, data):
        """
        This raises an exception of the appropriate type if the set of
        validation mechanisms are not met, and otherwise returns True.
        """
        for validator in self.validators:
            validator(data)
        # If we don't get an exception, we're good to go
        return True

    def get_dependencies(self, *args, **kwargs):
        """
        This returns a list of names of fields that this field depends on.
        """
        e = FieldDetector(*args, **kwargs)
        if self._function.__name__ == '<lambda>':
            e.requested.append(self.name)
        else:
            e[self.name]
        return e

    _unit_registry = None
    @contextlib.contextmanager
    def unit_registry(self, data):
        old_registry = self._unit_registry
        if hasattr(data, 'unit_registry'):
            ur = data.unit_registry
        elif hasattr(data, 'ds'):
            ur = data.ds.unit_registry
        else:
            ur = None
        self._unit_registry = ur
        yield
        self._unit_registry = old_registry

    def __call__(self, data):
        """ Return the value of the field in a given *data* object. """
        self.check_available(data)
        original_fields = data.keys() # Copy
        if self._function is NullFunc:
            raise RuntimeError(
                "Something has gone terribly wrong, _function is NullFunc " +
                "for %s" % (self.name,))
        with self.unit_registry(data):
            dd = self._function(self, data)
        for field_name in data.keys():
            if field_name not in original_fields:
                del data[field_name]
        return dd

    def get_source(self):
        """
        Return a string containing the source of the function (if possible.)
        """
        return inspect.getsource(self._function)

    def get_label(self, projected=False):
        """
        Return a data label for the given field, including units.
        """
        name = self.name[1]
        if self.display_name is not None:
            name = self.display_name

        # Start with the field name
        data_label = r"$\rm{%s}" % name

        # Grab the correct units
        if projected:
            raise NotImplementedError
        else:
            if self.ds is not None:
                units = Unit(self.units, registry=self.ds.unit_registry)
            else:
                units = Unit(self.units)
        # Add unit label
        if not units.is_dimensionless:
            data_label += r"\ \ (%s)" % (units.latex_representation())

        data_label += r"$"
        return data_label

    def __repr__(self):
        if PY2:
            func_name = self._function.func_name
        else:
            func_name = self._function.__name__
        if self._function == NullFunc:
            s = "On-Disk Field "
        elif func_name == "_TranslationFunc":
            s = "Alias Field for \"%s\" " % (self._function.alias_name,)
        else:
            s = "Derived Field "
        if isinstance(self.name, tuple):
            s += "(%s, %s): " % self.name
        else:
            s += "%s: " % (self.name)
        s += "(units: %s" % self.units
        if self.display_name is not None:
            s += ", display_name: '%s'" % (self.display_name)
        if self.particle_type:
            s += ", particle field"
        s += ")"
        return s

class FieldValidator(object):
    pass

class ValidateParameter(FieldValidator):
    def __init__(self, parameters, parameter_values=None):
        """
        This validator ensures that the dataset has a given parameter.

        If *parameter_values* is supplied, this will also ensure that the field
        is available for all permutations of the field parameter.
        """
        FieldValidator.__init__(self)
        self.parameters = ensure_list(parameters)
        self.parameter_values = parameter_values
    def __call__(self, data):
        doesnt_have = []
        if self.parameter_values is not None:
            if isinstance(data, FieldDetector):
                raise NeedsParameterValue(self.parameter_values)
        for p in self.parameters:
            if not data.has_field_parameter(p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsParameter(doesnt_have)
        return True

class ValidateDataField(FieldValidator):
    def __init__(self, field):
        """
        This validator ensures that the output file has a given data field stored
        in it.
        """
        FieldValidator.__init__(self)
        self.fields = ensure_list(field)
    def __call__(self, data):
        doesnt_have = []
        if isinstance(data, FieldDetector): return True
        for f in self.fields:
            if f not in data.index.field_list:
                doesnt_have.append(f)
        if len(doesnt_have) > 0:
            raise NeedsDataField(doesnt_have)
        return True

class ValidateProperty(FieldValidator):
    def __init__(self, prop):
        """
        This validator ensures that the data object has a given python attribute.
        """
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
        """
        This validator ensures that the data handed to the field is of spatial
        nature -- that is to say, 3-D.
        """
        FieldValidator.__init__(self)
        self.ghost_zones = ghost_zones
        self.fields = fields
    def __call__(self, data):
        # When we say spatial information, we really mean
        # that it has a three-dimensional data structure
        #if isinstance(data, FieldDetector): return True
        if not getattr(data, '_spatial', False):
            raise NeedsGridType(self.ghost_zones,self.fields)
        if self.ghost_zones <= data._num_ghost_zones:
            return True
        raise NeedsGridType(self.ghost_zones,self.fields)

class ValidateGridType(FieldValidator):
    def __init__(self):
        """
        This validator ensures that the data handed to the field is an actual
        grid patch, not a covering grid of any kind.
        """
        FieldValidator.__init__(self)
    def __call__(self, data):
        # We need to make sure that it's an actual AMR grid
        if isinstance(data, FieldDetector): return True
        if getattr(data, "_type_name", None) == 'grid': return True
        raise NeedsOriginalGrid()

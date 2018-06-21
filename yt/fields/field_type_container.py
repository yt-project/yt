"""
A proxy object for field descriptors, usually living as ds.fields.
"""

import weakref
from yt.extern.six import add_metaclass, string_types
from yt.fields.derived_field import \
    DerivedField

class FieldTypeContainer(object):
    def __init__(self, ds):
        self.ds = weakref.proxy(ds)

    def __getattr__(self, attr):
        ds = self.__getattribute__('ds')
        fnc = FieldNameContainer(ds, attr)
        if len(dir(fnc)) == 0:
            return self.__getattribute__(attr)
        return fnc

    _field_types = None
    @property
    def field_types(self):
        if self._field_types is None:
            self._field_types = set(t for t, n in self.ds.field_info)
        return self._field_types

    def __dir__(self):
        return list(self.field_types)

    def __iter__(self):
        for ft in self.field_types:
            fnc = FieldNameContainer(self.ds, ft)
            if len(dir(fnc)) == 0:
                yield self.__getattribute__(ft)
            else:
                yield fnc

    def __contains__(self, obj):
        ob = None
        if isinstance(obj, FieldNameContainer):
            ob = obj.field_type
        elif isinstance(obj, string_types):
            ob = obj

        return ob in self.field_types

class FieldNameContainer(object):
    def __init__(self, ds, field_type):
        self.ds = ds
        self.field_type = field_type

    def __getattr__(self, attr):
        ft = self.__getattribute__("field_type")
        ds = self.__getattribute__("ds")
        if (ft, attr) not in ds.field_info:
            return self.__getattribute__(attr)
        return ds.field_info[ft, attr]

    def __dir__(self):
        return [n for t, n in self.ds.field_info
                if t == self.field_type]

    def __iter__(self):
        for t, n in self.ds.field_info:
            if t == self.field_type:
                yield self.ds.field_info[t, n]

    def __contains__(self, obj):
        if isinstance(obj, DerivedField):
            if self.field_type == obj.name[0] and obj.name in self.ds.field_info:
                # e.g. from a completely different dataset
                if self.ds.field_info[obj.name] is not obj:
                    return False
                return True
        elif isinstance(obj, tuple):
            if self.field_type == obj[0] and obj in self.ds.field_info:
                return True
        elif isinstance(obj, string_types):
            if (self.field_type, obj) in self.ds.field_info:
                return True
        return False



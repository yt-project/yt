"""
A base-class representing an astrophysical object

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

astro_object_registry = {}

class AstrophysicalObject(object):
    # No _type_name
    _skip_add = False

    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_type_name") and not cls._skip_add:
                astro_object_registry[cls._type_name] = cls
            cls.identification_methods = {}
            cls.correlation_methods = {}

    def _lookup_object(self, obj_name):
        if obj_name not in astro_object_registry:
            raise KeyError(obj_name)
        return astro_object_registry[obj_name]

    def correlate(self, other_collection, correlation_name):
        pass

    def __init__(self, data_source):
        self.objects = {}
        # We mandate that every object have a corresponding YTSelectionContainer3D source
        # affiliated with it.
        self.data_source = data_source

    def find(self, obj_name, identification_name, *args, **kwargs):
        obj = self._lookup_object(obj_name)
        if callable(identification_name):
            identification_method = identification_name
        else:
            if identification_name not in obj.identification_methods:
                raise KeyError(identification_name)
            identification_method = \
                obj.identification_methods[identification_name]
        new_objs = identification_method(self, *args, **kwargs)
        setattr(self, obj_name, new_objs)
        self.objects[obj_name] = new_objs
        return new_objs

    def correlate(self, other_set, correlation_name, *args, **kwargs):
        if callable(correlation_name):
            correlation_method = correlation_name
        else:
            if correlation_name not in self.correlation_methods:
                raise KeyError(correlation_name)
            correlation_method = self.correlation_methods[correlation_name]
        linked_objs = correlation_method(self, *args, **kwargs)
        return linked_objs

def correlation_method(obj_name, link_name):
    def passthrough(func):
        obj = astro_object_registry[obj_name]
        obj.correlation_methods[link_name] = func
        return func
    return passthrough

def identification_method(obj_name, id_name):
    def passthrough(func):
        obj = astro_object_registry[obj_name]
        obj.identification_methods[id_name] = func
        return func
    return passthrough

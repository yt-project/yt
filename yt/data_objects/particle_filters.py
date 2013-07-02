"""
This is a library for defining and using particle filters.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

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

import numpy as np
import copy
from contextlib import contextmanager

from yt.data_objects.field_info_container import \
    NullFunc, TranslationFunc
from yt.utilities.exceptions import YTIllDefinedFilter
from yt.funcs import *

# One to many mapping
filter_registry = defaultdict(list)

class DummyFieldInfo(object):
    particle_type = True
dfi = DummyFieldInfo()

class ParticleFilter(object):
    def __init__(self, name, function, requires, filtered_type):
        self.name = name
        self.function = function
        self.requires = requires[:]
        self.filtered_type = filtered_type

    @contextmanager
    def apply(self, dobj):
        with dobj._chunked_read(dobj._current_chunk):
            with dobj._field_type_state(self.filtered_type, dfi):
                # We won't be storing the field data from the whole read, so we
                # start by filtering now.
                filter = self.function(self, dobj)
                yield
                # Retain a reference here, and we'll filter all appropriate fields
                # later.
                fd = dobj.field_data
        for f, tr in fd.items():
            if f[0] != self.filtered_type: continue
            if tr.shape != filter.shape:
                raise YTIllDefinedFilter(self, tr.shape, filter.shape)
            dobj.field_data[self.name, f[1]] = tr[filter]

    def available(self, field_list):
        # Note that this assumes that all the fields in field_list have the
        # same form as the 'requires' attributes.  This won't be true if the
        # fields are implicitly "all" or something.
        return all((self.filtered_type, field) in field_list for field in self.requires)

    def wrap_func(self, field_name, old_fi):
        new_fi = copy.copy(old_fi)
        new_fi.name = (self.filtered_type, field_name[1])
        return new_fi

def add_particle_filter(name, function, requires = None, filtered_type = "all"):
    if requires is None: requires = []
    filter = ParticleFilter(name, function, requires, filtered_type)
    filter_registry[name].append(filter)

def particle_filter(name, requires = None, filtered_type = "all"):
    def _pfilter(func):
        add_particle_filter(name, func, requires, filtered_type)
        return func
    return _pfilter

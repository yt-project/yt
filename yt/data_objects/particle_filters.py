"""
This is a library for defining and using particle filters.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import copy
from contextlib import contextmanager

from yt.fields.field_info_container import \
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
            if tr.shape != filter.shape and tr.shape[0] != filter.shape[0]:
                raise YTIllDefinedFilter(self, tr.shape, filter.shape)
            elif filter.size == 0:
                # Filtering empty set.  This keeps our dimensions correct.
                # Otherwise we end up with out-of-axis and shape problems.
                d = tr.copy() 
            elif len(tr.shape) > len(filter.shape):
                # Filter must always be 1D
                d = tr[filter,:]
            else:
                d = tr[filter]
            dobj.field_data[self.name, f[1]] = d

    def available(self, field_list):
        # Note that this assumes that all the fields in field_list have the
        # same form as the 'requires' attributes.  This won't be true if the
        # fields are implicitly "all" or something.
        return all((self.filtered_type, field) in field_list for field in self.requires)

    def wrap_func(self, field_name, old_fi):
        new_fi = copy.copy(old_fi)
        new_fi.name = (self.filtered_type, field_name[1])
        if old_fi._function == NullFunc:
            new_fi._function = TranslationFunc(old_fi.name)
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

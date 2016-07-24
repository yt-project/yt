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

import copy
from collections import defaultdict

from contextlib import contextmanager

from yt.fields.field_info_container import \
    NullFunc, TranslationFunc
from yt.utilities.exceptions import YTIllDefinedFilter

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


def add_particle_filter(name, function, requires=None, filtered_type="all"):
    r"""Create a new particle filter in the global namespace of filters

    A particle filter is a short name that corresponds to an algorithm for
    filtering a set of particles into a subset.  This is useful for creating new
    particle types based on a cut on a particle field, such as particle mass, ID
    or type.

    .. note::
       Alternatively, you can make use of the
       :func:`~yt.data_objects.particle_filters.particle_filter` decorator to
       define a new particle filter.

    Parameters
    ----------
    name : string
        The name of the particle filter.  New particle fields with particle type
        set by this name will be added to any dataset that enables this particle
        filter.
    function : reference to a function
        The function that defines the particle filter.  The function should
        accept two arguments: a reference to a particle filter object and a
        reference to an abstract yt data object.  See the example below.
    requires : a list of field names
        A list of field names required by the particle filter definition.
    filtered_type : string
        The name of the particle type to be filtered.

    Example
    -------

    >>> import yt

    >>> def _stars(pfilter, data):
    ...     return data[(pfilter.filtered_type, 'particle_type')] == 2

    >>> yt.add_particle_filter("stars", function=_stars, filtered_type='all',
    ...                        requires=["particle_type"])

    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> ds.add_particle_filter('stars')
    >>> ad = ds.all_data()
    >>> print (ad['stars', 'particle_mass'])
    [  1.68243760e+38   1.65690882e+38   1.65813321e+38 ...,   2.04238266e+38
       2.04523901e+38   2.04770938e+38] g

    """
    if requires is None:
        requires = []
    filter = ParticleFilter(name, function, requires, filtered_type)
    filter_registry[name].append(filter)


def particle_filter(name=None, requires=None, filtered_type='all'):
    r"""A decorator that adds a new particle filter

    A particle filter is a short name that corresponds to an algorithm for
    filtering a set of particles into a subset.  This is useful for creating new
    particle types based on a cut on a particle field, such as particle mass, ID
    or type.

    .. note::
       Alternatively, you can make use of the
       :func:`~yt.data_objects.particle_filters.add_particle_filter` function
       to define a new particle filter using a more declarative syntax.

    Parameters
    ----------
    name : string
        The name of the particle filter.  New particle fields with particle type
        set by this name will be added to any dataset that enables this particle
        filter.  If not set, the name will be inferred from the name of the
        filter function.
    function : reference to a function
        The function that defines the particle filter.  The function should
        accept two arguments: a reference to a particle filter object and a
        reference to an abstract yt data object.  See the example below.
    requires : a list of field names
        A list of field names required by the particle filter definition.
    filtered_type : string
        The name of the particle type to be filtered.

    Example
    -------

    >>> import yt

    >>> # define a filter named "stars"
    >>> @yt.particle_filter(requires=["particle_type"], filtered_type='all')
    >>> def stars(pfilter, data):
    ...     return data[(pfilter.filtered_type, 'particle_type')] == 2

    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> ds.add_particle_filter('stars')
    >>> ad = ds.all_data()
    >>> print (ad['stars', 'particle_mass'])
    [  1.68243760e+38   1.65690882e+38   1.65813321e+38 ...,   2.04238266e+38
       2.04523901e+38   2.04770938e+38] g

    """
    def wrapper(function):
        if name is None:
            used_name = function.__name__
        else:
            used_name = name
        return add_particle_filter(used_name, function, requires, filtered_type)
    return wrapper

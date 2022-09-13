import copy
from contextlib import contextmanager
from typing import Dict

from yt.fields.field_info_container import NullFunc, TranslationFunc
from yt.funcs import mylog
from yt.utilities.exceptions import YTIllDefinedFilter

# One to one mapping
filter_registry: Dict[str, "ParticleFilter"] = {}


class DummyFieldInfo:
    particle_type = True
    sampling_type = "particle"


dfi = DummyFieldInfo()


class ParticleFilter:
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
            if f[0] != self.filtered_type:
                continue
            if tr.shape != filter.shape and tr.shape[0] != filter.shape[0]:
                raise YTIllDefinedFilter(self, tr.shape, filter.shape)
            else:
                d = tr[filter]
            dobj.field_data[self.name, f[1]] = d

    def available(self, field_list):
        # Note that this assumes that all the fields in field_list have the
        # same form as the 'requires' attributes.  This won't be true if the
        # fields are implicitly "all" or something.
        return all((self.filtered_type, field) in field_list for field in self.requires)

    def missing(self, field_list):
        return list(
            (self.filtered_type, field)
            for field in self.requires
            if (self.filtered_type, field) not in field_list
        )

    def wrap_func(self, field_name, old_fi):
        new_fi = copy.copy(old_fi)
        new_fi.name = (self.name, field_name[1])
        if old_fi._function == NullFunc:
            new_fi._function = TranslationFunc(old_fi.name)
        # Marking the field as inherited
        new_fi._inherited_particle_filter = True
        return new_fi


def add_particle_filter(name, function, requires=None, filtered_type="all"):
    r"""Create a new particle filter in the global namespace of filters

    A particle filter is a short name that corresponds to an algorithm for
    filtering a set of particles into a subset.  This is useful for creating new
    particle types based on a cut on a particle field, such as particle mass, ID
    or type. After defining a new filter, it still needs to be added to the
    dataset by calling
    :func:`~yt.data_objects.static_output.add_particle_filter`.

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

    Examples
    --------
    >>> import yt

    >>> def _stars(pfilter, data):
    ...     return data[(pfilter.filtered_type, "particle_type")] == 2

    >>> yt.add_particle_filter(
    ...     "stars", function=_stars, filtered_type="all", requires=["particle_type"]
    ... )

    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ds.add_particle_filter("stars")
    >>> ad = ds.all_data()
    >>> print(ad["stars", "particle_mass"])
    [  1.68243760e+38   1.65690882e+38   1.65813321e+38 ...,   2.04238266e+38
       2.04523901e+38   2.04770938e+38] g

    """
    if requires is None:
        requires = []
    filter = ParticleFilter(name, function, requires, filtered_type)
    if filter_registry.get(name, None) is not None:
        mylog.warning("The %s particle filter already exists. Overriding.", name)
    filter_registry[name] = filter


def particle_filter(name=None, requires=None, filtered_type="all"):
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
    requires : a list of field names
        A list of field names required by the particle filter definition.
    filtered_type : string
        The name of the particle type to be filtered.

    Examples
    --------
    >>> import yt

    >>> # define a filter named "stars"
    >>> @yt.particle_filter(requires=["particle_type"], filtered_type="all")
    ... def stars(pfilter, data):
    ...     return data[(pfilter.filtered_type, "particle_type")] == 2

    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ds.add_particle_filter("stars")
    >>> ad = ds.all_data()
    >>> print(ad["stars", "particle_mass"])
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

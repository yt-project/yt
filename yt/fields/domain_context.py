import abc
from typing import Tuple

from yt._typing import FieldKey

domain_context_registry = {}


class DomainContext(abc.ABC):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            domain_context_registry[name] = cls

    _known_fluid_fields: Tuple[FieldKey, ...]
    _known_particle_fields: Tuple[FieldKey, ...]

    def __init__(self, ds):
        self.ds = ds

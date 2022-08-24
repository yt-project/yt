import abc
from typing import Tuple

domain_context_registry = {}


class DomainContext(abc.ABC):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            domain_context_registry[name] = cls

    _known_fluid_fields: Tuple[Tuple[str, str], ...]
    _known_particle_fields: Tuple[Tuple[str, str], ...]

    def __init__(self, ds):
        self.ds = ds

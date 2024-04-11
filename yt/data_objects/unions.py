from abc import ABC, abstractmethod

from more_itertools import always_iterable


class Union(ABC):
    @property
    @abstractmethod
    def _union_type(self) -> str: ...

    def __init__(self, name, sub_types):
        self.name = name
        self.sub_types = list(always_iterable(sub_types))

    def __iter__(self):
        yield from self.sub_types

    def __repr__(self):
        return f"{self._union_type.capitalize()} Union: '{self.name}' composed of: {self.sub_types}"


class MeshUnion(Union):
    _union_type = "mesh"


class ParticleUnion(Union):
    _union_type = "particle"

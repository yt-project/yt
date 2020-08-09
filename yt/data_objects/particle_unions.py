from typing import List

from .unions import Union


class ParticleUnion(Union):
    _union_type = "particle"

    def __init__(self, name: str, sub_types: List[str]) -> None:
        super(ParticleUnion, self).__init__(name, sub_types)

from .unions import Union


class ParticleUnion(Union):
    _union_type = "particle"

    def __init__(self, name, sub_types):
        super().__init__(name, sub_types)

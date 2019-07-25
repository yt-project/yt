from .unions import Union

class ParticleUnion(Union):
    def __init__(self, name, sub_types):
        super(ParticleUnion, self).__init__(name, sub_types)

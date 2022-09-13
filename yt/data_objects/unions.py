from more_itertools import always_iterable


class Union:
    _union_type = ""

    def __init__(self, name, sub_types):
        self.name = name
        self.sub_types = list(always_iterable(sub_types))

    def __iter__(self):
        yield from self.sub_types

    def __repr__(self):
        return "{} Union: '{}' composed of: {}".format(
            self._union_type.capitalize(), self.name, self.sub_types
        )


class MeshUnion(Union):
    _union_type = "mesh"

    def __init__(self, name, sub_types):
        super().__init__(name, sub_types)

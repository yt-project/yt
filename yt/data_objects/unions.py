from yt.funcs import ensure_list

class Union(object):
    _union_type = ""
    def __init__(self, name, sub_types):
        self.name = name
        self.sub_types = ensure_list(sub_types)

    def __iter__(self):
        for st in self.sub_types:
            yield st

    def __repr__(self):
        return "{} Union: '{}' composed of: {}".format(
            self._union_type.capitalize(), self.name, self.sub_types)

class MeshUnion(Union):
    _union_type = "mesh"
    def __init__(self, name, sub_types):
        super(MeshUnion, self).__init__(name, sub_types)

from yt.funcs import ensure_list

class Union(object):
    def __init__(self, name, sub_types):
        self.name = name
        self.sub_types = ensure_list(sub_types)

    def __iter__(self):
        for st in self.sub_types:
            yield st

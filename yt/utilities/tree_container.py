class TreeContainer:
    r"""A recursive data container for things like merger trees and
    clump-finder trees.

    """
    _child_attr = "children"

    def __init__(self):
        setattr(self, self._child_attr, None)

    def __iter__(self):
        yield self
        children = getattr(self, self._child_attr)
        if children is None:
            return
        for child in children:
            yield from child

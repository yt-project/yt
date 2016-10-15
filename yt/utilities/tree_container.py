"""
TreeContainer class and member functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

class TreeContainer(object):
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
            for a_node in child:
                yield a_node

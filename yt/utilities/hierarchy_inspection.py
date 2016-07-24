# -*- coding: utf-8 -*-
import inspect
from collections import Counter
from yt.extern.six.moves import reduce


def find_lowest_subclasses(candidates):
    """
    This function takes a list of classes, and returns only the ones that are
    are not super classes of any others in the list. i.e. the ones that are at
    the bottom of the specified mro of classes.

    Parameters
    ----------
    candidates : iterable
        An interable object that is a collection of classes to find the lowest
        subclass of.

    Returns
    -------
    result : list
        A list of classes which are not super classes for any others in
        candidates.
    """

    # If there is only one input, the input candidate is always the
    # lowest class
    if len(candidates) == 1:
        return candidates
    elif len(candidates) == 0:
        return []

    mros = [inspect.getmro(c) for c in candidates]

    counters = [Counter(mro) for mro in mros]

    if len(counters) == 0:
        return counters

    count = reduce(lambda x, y: x + y, counters)

    return [x for x in count.keys() if count[x] == 1]

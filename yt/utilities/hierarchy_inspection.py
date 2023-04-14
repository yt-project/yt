import inspect
from collections import Counter
from typing import List, Type

from more_itertools import flatten


def find_lowest_subclasses(candidates: List[Type]) -> List[Type]:
    """
    This function takes a list of classes, and returns only the ones that are
    are not super classes of any others in the list. i.e. the ones that are at
    the bottom of the specified mro of classes.

    Parameters
    ----------
    candidates : Iterable
        An iterable object that is a collection of classes to find the lowest
        subclass of.

    Returns
    -------
    result : list
        A list of classes which are not super classes for any others in
        candidates.
    """
    count = Counter(flatten(inspect.getmro(c) for c in candidates))
    return [x for x in candidates if count[x] == 1]

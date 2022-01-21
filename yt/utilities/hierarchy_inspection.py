import inspect
from collections import Counter
from functools import reduce
from typing import List, Optional, Type


def find_lowest_subclasses(
    candidates: List[Type], *, hint: Optional[str] = None
) -> List[Type]:
    """
    This function takes a list of classes, and returns only the ones that are
    are not super classes of any others in the list. i.e. the ones that are at
    the bottom of the specified mro of classes.

    Parameters
    ----------
    candidates : Iterable
        An iterable object that is a collection of classes to find the lowest
        subclass of.

    hint : str, optional
        Only keep candidates classes that have `hint` in their name (case insensitive)

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
        return []

    count = reduce(lambda x, y: x + y, counters)

    retv = [x for x in count.keys() if count[x] == 1]
    if hint is not None:
        retv = [x for x in retv if hint.lower() in x.__name__.lower()]
    return retv

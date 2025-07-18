import inspect
from collections import Counter
from typing import TypeVar

from more_itertools import flatten

from yt.data_objects.static_output import Dataset
from yt.utilities.object_registries import output_type_registry

T = TypeVar("T")


def find_lowest_subclasses(candidates: list[type[T]]) -> list[type[T]]:
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


def get_classes_with_missing_requirements() -> dict[type[Dataset], list[str]]:
    # We need a function here rather than an global constant registry because:
    # - computation should be delayed until needed so that the result is independent of import order
    # - tests might (temporarily) mutate output_type_registry
    return {
        cls: missing
        for cls in sorted(output_type_registry.values(), key=lambda c: c.__name__)
        if (missing := cls._missing_load_requirements())
    }

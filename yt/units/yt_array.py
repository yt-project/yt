from unyt.array import *

from yt.funcs import array_like_field
from yt.units import YTArray, YTQuantity


def display_ytarray(arr):
    r"""
    Display a YTArray in a Jupyter widget that enables unit switching.

    The array returned by this function is read-only, and only works with
    arrays of size 3 or lower.

    Parameters
    ----------
    arr : YTArray
        The Array to display; must be of size 3 or lower.

    Examples
    --------
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> display_ytarray(ds.domain_width)
    """
    if arr.size > 3:
        raise YTArrayTooLargeToDisplay(arr.size, 3)
    import ipywidgets

    unit_registry = arr.units.registry
    equiv = unit_registry.list_same_dimensions(arr.units)
    dropdown = ipywidgets.Dropdown(options=sorted(equiv), value=str(arr.units))

    def arr_updater(arr, texts):
        def _value_updater(change):
            arr2 = arr.in_units(change["new"])
            if arr2.shape == ():
                arr2 = [arr2]
            for v, t in zip(arr2, texts):
                t.value = str(v.value)

        return _value_updater

    if arr.shape == ():
        arr_iter = [arr]
    else:
        arr_iter = arr
    texts = [ipywidgets.Text(value=str(_.value), disabled=True) for _ in arr_iter]
    dropdown.observe(arr_updater(arr, texts), names="value")
    return ipywidgets.HBox(texts + [dropdown])


def _wrap_display_ytarray(arr):
    from IPython.core.display import display

    display(display_ytarray(arr))

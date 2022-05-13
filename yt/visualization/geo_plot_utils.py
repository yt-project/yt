from types import FunctionType
from typing import Any, Dict, Optional, Tuple

valid_transforms: Dict[str, FunctionType] = {}

transform_list = [
    "PlateCarree",
    "LambertConformal",
    "LambertCylindrical",
    "Mercator",
    "Miller",
    "Mollweide",
    "Orthographic",
    "Robinson",
    "Stereographic",
    "TransverseMercator",
    "InterruptedGoodeHomolosine",
    "RotatedPole",
    "OSGB",
    "EuroPP",
    "Geostationary",
    "Gnomonic",
    "NorthPolarStereo",
    "OSNI",
    "SouthPolarStereo",
    "AlbersEqualArea",
    "AzimuthalEquidistant",
    "Sinusoidal",
    "UTM",
    "NearsidePerspective",
    "LambertAzimuthalEqualArea",
]


def cartopy_importer(transform_name):
    r"""Convenience function to import cartopy projection types"""

    def _func(*args, **kwargs):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        return getattr(cartopy.crs, transform_name)(*args, **kwargs)

    return _func


def get_mpl_transform(mpl_proj) -> Optional[FunctionType]:
    r"""This returns an instantiated transform function given a transform
    function name and arguments.

    Parameters
    ----------
    mpl_proj : string or tuple
        the matplotlib projection type. Can take the form of string or tuple.

    Examples
    --------

    >>> get_mpl_transform("PlateCarree")

    >>> get_mpl_transform(
    ...     ("RotatedPole", (), {"pole_latitude": 37.5, "pole_longitude": 177.5})
    ... )

    """
    # first check to see if the transform dict is empty, if it is fill it with
    # the cartopy functions
    if not valid_transforms:
        for mpl_transform in transform_list:
            valid_transforms[mpl_transform] = cartopy_importer(mpl_transform)

    # check to see if mpl_proj is a string or tuple, and construct args and
    # kwargs to pass to cartopy function based on that.
    key: Optional[str] = None
    args: Tuple = ()
    kwargs: Dict[str, Any] = {}
    if isinstance(mpl_proj, str):
        key = mpl_proj
        instantiated_func = valid_transforms[key](*args, **kwargs)
    elif isinstance(mpl_proj, tuple):
        if len(mpl_proj) == 2:
            key, args = mpl_proj
            kwargs = {}
        elif len(mpl_proj) == 3:
            key, args, kwargs = mpl_proj
        else:
            raise ValueError(f"Expected a tuple with len 2 or 3, received {mpl_proj}")
        if not isinstance(key, str):
            raise TypeError(
                f"Expected a string a the first element in mpl_proj, got {key!r}"
            )
        instantiated_func = valid_transforms[key](*args, **kwargs)
    elif hasattr(mpl_proj, "globe"):
        # cartopy transforms have a globe method associated with them
        key = mpl_proj
        instantiated_func = mpl_proj
    elif hasattr(mpl_proj, "set_transform"):
        # mpl axes objects have a set_transform method, so we'll check if that
        # exists on something passed in.
        key = mpl_proj
        instantiated_func = mpl_proj
    elif hasattr(mpl_proj, "name"):
        # last we'll check if the transform is in the list of registered
        # projections in matplotlib.
        from matplotlib.projections import get_projection_names

        registered_projections = get_projection_names()
        if mpl_proj.name in registered_projections:
            key = mpl_proj
            instantiated_func = mpl_proj
        else:
            key = None

    # build in a check that if none of the above options are satisfied by what
    # the user passes that None is returned for the instantiated function
    if key is None:
        instantiated_func = None
    return instantiated_func

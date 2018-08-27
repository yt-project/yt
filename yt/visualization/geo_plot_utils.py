valid_transforms={}

transform_list = ['PlateCarree', 'LambertConformal', 'LambertCylindrical',
                  'Mercator', 'Miller', 'Mollweide', 'Orthographic',
                  'Robinson', 'Stereographic', 'TransverseMercator',
                  'InterruptedGoodeHomolosine', 'RotatedPole', 'OSGB',
                  'EuroPP', 'Geostationary', 'Gnomonic', 'NorthPolarStereo',
                  'OSNI', 'SouthPolarStereo', 'AlbersEqualArea',
                  'AzimuthalEquidistant', 'Sinusoidal', 'UTM',
                  'NearsidePerspective', 'LambertAzimuthalEqualArea']

def cartopy_importer(transform_name):
    r"""Convenience function to import cartopy projection types
    """

    def _func(*args, **kwargs):
        from yt.utilities.on_demand_imports import _cartopy as cartopy
        return getattr(cartopy.crs, transform_name)(*args, **kwargs)
    return _func

def get_mpl_transform(mpl_proj):
    r"""This returns an instantiated transform function given a transform
    function name and arguments.

    Parameters
    ----------
    mpl_proj : string or tuple
        the matplotlib projection type. Can take the form of string or tuple.

    Examples
    --------

    >>> get_mpl_transform("PlateCarree")

    >>> get_mpl_transform(('RotatedPole', (), {'pole_latitude':37.5,
    ...                    'pole_longitude':177.5}))

    """
    # first check to see if the tranform dict is empty, if it is fill it with
    # the cartopy functions
    if not valid_transforms:
        for mpl_transform in transform_list:
            valid_transforms[mpl_transform] = cartopy_importer(mpl_transform)

    # check to see if mpl_proj is a string or tuple, and construct args and
    # kwargs to pass to cartopy function based on that.
    key = None
    if isinstance (mpl_proj, str):
        key = mpl_proj
        args = ()
        kwargs = {}
        instantiated_func = valid_transforms[key](*args, **kwargs)
    elif isinstance (mpl_proj, tuple):
        if len(mpl_proj) == 2:
            key, args = mpl_proj
            kwargs = {}
        elif len(mpl_proj) == 3:
            key, args, kwargs = mpl_proj
        instantiated_func = valid_transforms[key](*args, **kwargs)
    if key is None:
        instantiated_func = None
    return instantiated_func

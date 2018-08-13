def cartopy_importer(transform_name):
        def _func(*args, **kwargs):
            import cartopy.crs
            return getattr(cartopy.crs, transform_name)(*args, **kwargs)
        return _func


valid_transforms={}

for mpl_transform in ['PlateCarree', 'Mollweide', 'Orthographic', 'Robinson']:
    valid_transforms[mpl_transform] = cartopy_importer(mpl_transform)

def get_mpl_transform(mpl_proj):
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

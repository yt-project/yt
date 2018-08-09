import matplotlib
matplotlib.rc('contour', negative_linestyle='solid')

from matplotlib.backends.backend_agg import \
    FigureCanvasAgg

from matplotlib.backends.backend_pdf import \
    FigureCanvasPdf

from matplotlib.backends.backend_ps import \
    FigureCanvasPS

valid_transform_list = ['PlateCaree', 'Mollweide', 'Orthographic', 'Robinson']

valid_transforms={}

def cartopy_importer(transform_name):
        def _func(*args, **kwargs):
            import cartopy.crs
            if transform_name in valid_transform_list:
                return getattr(cartopy.crs, transform_name)(*args, **kwargs)
        return _func

def get_transform(mpl_proj):
    key = None
    if isinstance (mpl_proj, str):
        key = mpl_proj
        args = ()
        kwargs = {}
    elif isinstance (mpl_proj, tuple):
        if len(mpl_proj) == 2:
            key, args = mpl_proj
            kwargs = {}
        elif len(mpl_proj) == 3:
            key args, kwargs = mpl_proj
    if key is None:
        raise RuntimeError

instantiated_func = valid_transforms[key](*args **kwargs)

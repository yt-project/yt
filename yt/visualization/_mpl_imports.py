import matplotlib
matplotlib.rc('contour', negative_linestyle='solid')

import matplotlib.image
import matplotlib.ticker
import matplotlib.axes
import matplotlib.figure
import matplotlib._image
import matplotlib.colors
import matplotlib.colorbar
import matplotlib.cm
import matplotlib.collections

from matplotlib.backends.backend_agg import \
    FigureCanvasAgg

from matplotlib.backends.backend_pdf import \
    FigureCanvasPdf

from matplotlib.backends.backend_ps import \
    FigureCanvasPS

# Now we provide some convenience functions to get information about plots.
# With Matplotlib 0.98.x, the 'transforms' branch broke backwards
# compatibility.  Despite that, the various packagers are plowing ahead with
# packaging 0.98.x with new distributions of python software.  So I guess
# we have to support it.

_compatibility_functions = ["mpl_get_bounds","mpl_notify"]

_mpl98_mpl_get_bounds = lambda bbox: bbox.bounds
_mpl9x_mpl_get_bounds = lambda bbox: bbox.get_bounds()
_mpl98_mpl_notify = lambda im,cb: cb.update_bruteforce(im)
_mpl9x_mpl_notify = lambda im,cb: cb.notify(im)

# This next function hurts, because it relies on the fact that we're
# only differentiating between 0.9[01] and 0.98. And if happens to be
# 1.0, or any version with only 3 values, this should catch it.

try:
    _mpl_version = float(matplotlib.__version__[:4])
except:
    _mpl_version = float(matplotlib.__version__[:3])

if _mpl_version < 0.98:
    _prefix = '_mpl9x'
else:
    _prefix = '_mpl98'

for fn in _compatibility_functions:
    exec("%s = %s_%s" % (fn, _prefix, fn))

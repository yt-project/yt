"""
This is an interface to U{MatPlotLib <http://matplotlib.sf.net>} to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.raven import *
import AMRPixelize

try:
    import matplotlib
except ImportError:
    mylog.warning("matplotlib failed to import; all AMRPlot commands will fail")

import matplotlib.image
import matplotlib.axes
import matplotlib._image
import matplotlib.colors

class _AMRImage(matplotlib.image.NonUniformImage):
    _buff = None
    def make_image(self, magnification=1.0):
        if self._A is None:
            raise RuntimeError('You must first set the image array')
        if self._buff != None:
            del self._buff

        x0, y0, v_width, v_height = self.axes.viewLim.get_bounds()
        l, b, width, height = self.axes.bbox.get_bounds()
        width *= magnification
        height *= magnification
        buff = AMRPixelize.Pixelize(self._Ax, self._Ay, 
                                    self._Adx, self._Ady,
                                    self._A, int(height), int(width),
                                    x0, x0+v_width, y0, y0+v_height,
                                    )
        self.norm.autoscale(buff)
        if self.next_clim[0] is not None: self.norm.vmin = self.next_clim[0]
        if self.next_clim[1] is not None: self.norm.vmax = self.next_clim[1]

        Ax = (self.cmap(self.norm(buff))*255).astype(nT.UInt8)
        im = matplotlib._image.frombyte(Ax, 1)

        self.next_clim = (None, None)

        bg = matplotlib.colors.colorConverter.to_rgba(self.axes.get_frame().get_facecolor(), 0)
        im.set_bg(*bg)
        self._buff = buff
        del buff
        return im

    def set_data(self, x, y, dx, dy, d):
        x = na.asarray(x).astype(nT.Float64)
        y = na.asarray(y).astype(nT.Float64)
        dx = na.asarray(dx).astype(nT.Float64)
        dy = na.asarray(dy).astype(nT.Float64)
        d = na.asarray(d)
        ind=na.argsort(dx)
        self._A = d[ind][::-1]
        self._Ax = x[ind][::-1]
        self._Ay = y[ind][::-1]
        self._Adx = dx[ind][::-1]
        self._Ady = dy[ind][::-1]
        self._imcache = None
        self.next_clim = (None, None)

    def set_next_clim(self, vmin, vmax):
        self.next_clim = (vmin, vmax)



def amrshow(self, x, y, dx, dy, A,
           cmap = None,
           norm = None,
           aspect=None,
           interpolation=None,
           alpha=1.0,
           vmin = None,
           vmax = None,
           origin=None,
           extent=None,
           shape=None,
           filternorm=1,
           filterrad=4.0,
           imlim=None,
           **kwargs):
    """

    IMSHOW(X, cmap=None, norm=None, aspect=None, interpolation=None,
           alpha=1.0, vmin=None, vmax=None, origin=None, extent=None)

    IMSHOW(X) - plot image X to current axes, resampling to scale to axes
                size (X may be numarray/Numeric array or PIL image)

    IMSHOW(X, **kwargs) - Use keyword args to control image scaling,
    colormapping etc. See below for details


    Display the image in X to current axes.  X may be a float array, a
    UInt8 array or a PIL image. If X is an array, X can have the following
    shapes:

        MxN    : luminance (grayscale, float array only)

        MxNx3  : RGB (float or UInt8 array)

        MxNx4  : RGBA (float or UInt8 array)

    The value for each component of MxNx3 and MxNx4 float arrays should be
    in the range 0.0 to 1.0; MxN float arrays may be normalised.

    A matplotlib.image.AxesImage instance is returned

    The following kwargs are allowed:

      * cmap is a cm colormap instance, eg cm.jet.  If None, default to rc
        image.cmap value (Ignored when X has RGB(A) information)

      * aspect is one of: auto, equal, or a number.  If None, default to rc
        image.aspect value

      * interpolation is one of:

        'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36',
        'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
        'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc',
        'lanczos', 'blackman'

        if interpolation is None, default to rc
        image.interpolation.  See also th the filternorm and
        filterrad parameters

      * norm is a matplotlib.colors.Normalize instance; default is
        normalization().  This scales luminance -> 0-1 (only used for an
        MxN float array).

      * vmin and vmax are used to scale a luminance image to 0-1.  If
        either is None, the min and max of the luminance values will be
        used.  Note if you pass a norm instance, the settings for vmin and
        vmax will be ignored.

      * alpha = 1.0 : the alpha blending value

      * origin is either upper or lower, which indicates where the [0,0]
        index of the array is in the upper left or lower left corner of
        the axes.  If None, default to rc image.origin

      * extent is a data xmin, xmax, ymin, ymax for making image plots
        registered with data plots.  Default is the image dimensions
        in pixels

      * shape is for raw buffer images

      * filternorm is a parameter for the antigrain image resize
        filter.  From the antigrain documentation, if normalize=1,
        the filter normalizes integer values and corrects the
        rounding errors. It doesn't do anything with the source
        floating point values, it corrects only integers according
        to the rule of 1.0 which means that any sum of pixel
        weights must be equal to 1.0.  So, the filter function
        must produce a graph of the proper shape.

     * filterrad: the filter radius for filters that have a radius
       parameter, ie when interpolation is one of: 'sinc',
       'lanczos' or 'blackman'

    Additional kwargs are matplotlib.artist properties
    """

    if not self._hold: self.cla()

    if norm is not None: assert(isinstance(norm, matplotlib.colors.Normalize))
    if cmap is not None: assert(isinstance(cmap, matplotlib.colors.Colormap))
    #if aspect is None: aspect = rcParams['image.aspect']
    if aspect is None: aspect = 1.0
    self.set_aspect(aspect)
    im = _AMRImage(self, cmap, norm, extent)

    im.set_data(x, y, dx, dy, A)
    #im.set_alpha(alpha)
    self._set_artist_props(im)
    #if norm is None and shape is None:
    #    im.set_clim(vmin, vmax)
    if vmin is not None or vmax is not None:
        im.set_clim(vmin, vmax)
    else:
        im.autoscale()

    xmin, xmax, ymin, ymax = im.get_extent()

    corners = (xmin, ymin), (xmax, ymax)
    self.update_datalim(corners)
    if self._autoscaleon:
        self.set_xlim((xmin, xmax))
        self.set_ylim((ymin, ymax))
    self.images.append(im)

    return im

matplotlib.axes.Axes.amrshow = amrshow

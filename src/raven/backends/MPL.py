"""
This is an interface to U{MatPlotLib <http://matplotlib.sf.net>} to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.raven import *
from yt.funcs import *

# We only get imported if matplotlib was imported successfully

import _MPL

import matplotlib.image
import matplotlib.axes
import matplotlib._image
import matplotlib.colors
from matplotlib.backends.backend_agg import FigureCanvasAgg

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
        buff = _MPL.Pixelize(self._Ax, self._Ay, 
                                    self._Adx, self._Ady,
                                    self._Adata, int(height), int(width),
                                   (x0, x0+v_width, y0, y0+v_height),
                                    )
        self.norm.autoscale(buff)
        self.norm.vmin = buff.min()
        self.norm.vmax = buff.max()
        #print "NORM", self.norm.vmin, self.norm.vmax, buff.min(), buff.max()
        if self.next_clim[0] is not None: self.norm.vmin = self.next_clim[0]
        if self.next_clim[1] is not None: self.norm.vmax = self.next_clim[1]

        Ax = (self.cmap(self.norm(buff))*255).astype(nT.UInt8)
        im = matplotlib._image.frombyte(Ax, 1)

        self.next_clim = (None, None)

        bg = matplotlib.colors.colorConverter.to_rgba(self.axes.get_frame().get_facecolor(), 0)
        im.set_bg(*bg)
        self._buff = buff
        self._A = buff
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
        self._Adata = d[ind][::-1]
        self._Ax = x[ind][::-1]
        self._Ay = y[ind][::-1]
        self._Adx = dx[ind][::-1]
        self._Ady = dy[ind][::-1]
        self._imcache = None
        self.next_clim = (None, None)

    def set_next_clim(self, vmin, vmax):
        self.next_clim = (vmin, vmax)

    def setWidth(self, width, unit):
        pass



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

class LinkedAMRSubPlots:
    def __init__(self, linkWidth, linkZ, plots = []):
        pass
    def SaveFigure(self, filename, format):
        pass
    def setWidth(self, width, unit):
        pass
    def setZ(self, zmin, zmax):
        pass
    def __getitem__(self, item):
        pass

class ClassicThreePane:
    pass

def ClusterFilePlot(cls, x, y, xlog=None, ylog=None, fig=None, filename=None,
                    format="png", xbounds = None, ybounds = None):
    """
    
    """
    if not fig:
        fig = matplotlib.figure.Figure(figsize=(8,8))
        canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    if not iterable(cls):
        cls = [cls]
    if xlog == None:
        if lagos.CFfieldInfo.has_key(x):
            xlog = lagos.CFfieldInfo[x][2]
    if ylog == None:
        if lagos.CFfieldInfo.has_key(y):
            ylog = lagos.CFfieldInfo[y][2]
    if xlog and ylog:
        pp=ax.loglog
    elif xlog and not ylog:
        pp=ax.semilogx
    elif ylog and not xlog:
        pp=ax.semilogy
    else:
        pp=ax.plot

    fig.hold(True)
    for cl in cls:
        pp(cl[x],cl[y])
    if lagos.CFfieldInfo.has_key(x):
        ax.set_xlabel(lagos.CFfieldInfo[x][1])
    if lagos.CFfieldInfo.has_key(y):
        ax.set_ylabel(lagos.CFfieldInfo[y][1])
    if xbounds:
        ax.set_xlim(xbounds)
    if ybounds:
        ax.set_ylim(ybounds)
    if filename:
        canvas.print_figure(filename, format=format)
    else:
        return fig

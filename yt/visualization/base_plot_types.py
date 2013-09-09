"""
This is a place for base classes of the various plot types.


Authors:
 * Nathan Goldbaum 


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import matplotlib
import cStringIO
from ._mpl_imports import \
    FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS
from yt.funcs import \
    get_image_suffix, mylog


class PlotMPL(object):
    """A base class for all yt plots made using matplotlib.

    """
    def __init__(self, fsize, axrect, figure, axes):
        """Initialize PlotMPL class"""
        self._plot_valid = True
        if figure is None:
            self.figure = matplotlib.figure.Figure(figsize=fsize, frameon=True)
        else:
            figure.set_size_inches(fsize)
            self.figure = figure
        if axes is None:
            self.axes = self.figure.add_axes(axrect)
        else:
            axes.cla()
            self.axes = axes
        self.canvas = FigureCanvasAgg(self.figure)

    def save(self, name, mpl_kwargs, canvas=None):
        """Choose backend and save image to disk"""
        suffix = get_image_suffix(name)
        if suffix == '':
            suffix = '.png'
            name = "%s%s" % (name, suffix)

        mylog.info("Saving plot %s", name)

        if suffix == ".png":
            canvas = FigureCanvasAgg(self.figure)
        elif suffix == ".pdf":
            canvas = FigureCanvasPdf(self.figure)
        elif suffix in (".eps", ".ps"):
            canvas = FigureCanvasPS(self.figure)
        else:
            mylog.warning("Unknown suffix %s, defaulting to Agg", suffix)
            canvas = self.canvas

        canvas.print_figure(name, **mpl_kwargs)
        return name


class ImagePlotMPL(PlotMPL):
    """A base class for yt plots made using imshow

    """
    def __init__(self, fsize, axrect, caxrect, zlim, figure, axes, cax):
        """Initialize ImagePlotMPL class object"""
        PlotMPL.__init__(self, fsize, axrect, figure, axes)
        self.zmin, self.zmax = zlim
        if cax is None:
            self.cax = self.figure.add_axes(caxrect)
        else:
            cax.cla()
            cax.set_position(caxrect)
            self.cax = cax

    def _init_image(self, data, cbnorm, cmap, extent, aspect):
        """Store output of imshow in image variable"""
        if (cbnorm == 'log10'):
            norm = matplotlib.colors.LogNorm()
        elif (cbnorm == 'linear'):
            norm = matplotlib.colors.Normalize()
        self.image = self.axes.imshow(data, origin='lower', extent=extent,
                                      norm=norm, vmin=self.zmin, aspect=aspect,
                                      vmax=self.zmax, cmap=cmap)
        self.cb = self.figure.colorbar(self.image, self.cax)

    def _repr_png_(self):
        canvas = FigureCanvasAgg(self.figure)
        f = cStringIO.StringIO()
        canvas.print_figure(f)
        f.seek(0)
        return f.read()

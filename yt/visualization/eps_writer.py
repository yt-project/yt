"""
DualEPS: A class to combine bitmap compression and vector graphics



"""
from __future__ import absolute_import, print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import pyx
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from yt.config import \
    ytcfg
from yt.utilities.logger import ytLogger as mylog
from .plot_window import PlotWindow
from .profile_plotter import PhasePlot, ProfilePlot
from yt.units.yt_array import YTQuantity
from yt.units.unit_object import Unit

def convert_frac_to_tex(string):
    frac_pos = string.find(r'\frac')
    result = string[frac_pos+5:]
    level = [0]*len(result)
    clevel = 0
    for i in range(len(result)):
        if result[i] == "{":
            clevel += 1
        elif result[i] == "}":
            clevel -= 1
        level[i] = clevel
    div_pos = level.index(0)
    end_pos = level.index(0,div_pos+1)
    result = r'${' + result[:div_pos+1] + r'\over' + result[div_pos+1:end_pos] + \
             r'}$' + result[end_pos:]
    result = result.replace(r'\ ', r'\;')
    return result

def pyxize_label(string):
    frac_pos = string.find(r'\frac')
    if frac_pos >= 0:
        pre = string[:frac_pos]
        result = pre + convert_frac_to_tex(string)
    else:
        result = string
    result = result.replace('$', '')
    result = r'$' + result + r'$'
    return result


class DualEPS(object):
    def __init__(self, figsize=(12,12)):
        r"""Initializes the DualEPS class to which we can progressively add layers
        of vector graphics and compressed bitmaps.

        Parameters
        ----------
        figsize : tuple of floats
            The width and height of a single figure in centimeters.
        """
        pyx.unit.set(xscale=1.4)
        self.figsize = figsize
        self.canvas = None
        self.colormaps = None
        self.field = None
        self.axes_drawn = False

    def hello_world(self):
        r"""A simple test.
        """
        if self.canvas is None:
            self.canvas = pyx.canvas.canvas()
        p = pyx.path.line(0,0,1,1)
        self.canvas.stroke(p)
        self.canvas.text(0,0, "Hello world.")

#=============================================================================

    def return_field(self, plot):
        if isinstance(plot, (PlotWindow, PhasePlot)):
            return plot.plots.keys()[0]
        else:
            return None

#=============================================================================

    def axis_box(self, xrange=(0,1), yrange=(0,1), xlabel="", ylabel="",
                 xlog=False, ylog=False, xdata=None, ydata=None,
                 tickcolor=None, bare_axes=False,
                 pos=(0,0), xaxis_side=0, yaxis_side=0, size=None):
        r"""Draws an axis box in the figure.

        Parameters
        ----------
        xrange : tuple of floats
            The min and max of the x-axis
        yrange : tuple of floats
            The min and max of the y-axis
        xlabel : string
            Label for the x-axis
        ylabel : string
            Label for the y-axis
        xlog : boolean
            Flag to use a logarithmic x-axis
        ylog : boolean
            Flag to use a logarithmic y-axis
        tickcolor : `pyx.color.*.*`
            Color for the tickmarks.  Example: pyx.color.cmyk.black
        bare_axes : boolean
            Set to true to have no annotations or tick marks on all of the
            axes.
        pos : tuple of floats
            (x,y) position in centimeters of the origin in the figure
        xaxis_side : integer
            Set to 0 for the x-axis annotations to be on the left.  Set
            to 1 to print them on the right side.
        yaxis_side : integer
            Set to 0 for the y-axis annotations to be on the bottom.  Set
            to 1 to print them on the top.
        size : tuple of floats
            Size of axis box in units of figsize

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.save_fig()
        """
        
        if isinstance(xrange[0], YTQuantity):
            xrange = (xrange[0].value, xrange[1].value)
        if isinstance(yrange[0], YTQuantity):
            yrange = (yrange[0].value, yrange[1].value)
        if tickcolor is None:
            c1 = pyx.graph.axis.painter.regular\
                 (tickattrs=[pyx.color.cmyk.black])
            c2 = pyx.graph.axis.painter.regular\
                 (tickattrs=[pyx.color.cmyk.black], labelattrs=None)
        else:
            c1 = pyx.graph.axis.painter.regular(tickattrs=[tickcolor])
            c2 = pyx.graph.axis.painter.regular\
                 (tickattrs=[tickcolor], labelattrs=None)

        if size is None:
            psize = self.figsize
        else:
            psize = (size[0]*self.figsize[0], size[1]*self.figsize[1])

        xticklabels = True
        yticklabels = True
        if xaxis_side == 0:
            xleftlabel = xlabel
            xrightlabel = ""
            c1x = c1
            c2x = c2
        elif xaxis_side == 1:
            xleftlabel = ""
            xrightlabel = xlabel
            c1x = c2
            c2x = c1
        else:
            xticklabels = False
            xleftlabel = ""
            xrightlabel = ""
            c1x = c1
            c2x = c2
        if yaxis_side == 0:
            yleftlabel = ylabel
            yrightlabel = ""
            c1y = c1
            c2y = c2
        elif yaxis_side == 1:
            yleftlabel = ""
            yrightlabel = ylabel
            c1y = c2
            c2y = c1
        else:
            yticklabels = False
            yleftlabel = ""
            yrightlabel = ""
            c1y = c1
            c2y = c2
        

     
        if xlog:
            if xticklabels:
                xaxis = pyx.graph.axis.log(min=xrange[0],max=xrange[1],
                                           title=xleftlabel, painter=c1x)
                xaxis2 = pyx.graph.axis.log(min=xrange[0],max=xrange[1],
                                            title=xrightlabel, painter=c2x)
            else:
                xaxis = pyx.graph.axis.log(min=xrange[0],max=xrange[1],
                                           title=xleftlabel, painter=c1x,
                                           parter=None)
                xaxis2 = pyx.graph.axis.log(min=xrange[0],max=xrange[1],
                                            title=xrightlabel, painter=c2x,
                                            parter=None)
        else:
            if xticklabels:
                xaxis = pyx.graph.axis.lin(min=xrange[0],max=xrange[1],
                                           title=xleftlabel, painter=c1x)
                xaxis2 = pyx.graph.axis.lin(min=xrange[0],max=xrange[1],
                                            title=xrightlabel, painter=c2x)
            else:
                xaxis = pyx.graph.axis.lin(min=xrange[0],max=xrange[1],
                                           title=xleftlabel, painter=c1x,
                                           parter=None)
                xaxis2 = pyx.graph.axis.lin(min=xrange[0],max=xrange[1],
                                            title=xrightlabel, painter=c2x,
                                            parter=None)
        if ylog:
            if yticklabels:
                yaxis = pyx.graph.axis.log(min=yrange[0],max=yrange[1],
                                           title=yleftlabel, painter=c1y)
                yaxis2 = pyx.graph.axis.log(min=yrange[0],max=yrange[1],
                                            title=yrightlabel, painter=c2y)
            else:
                yaxis = pyx.graph.axis.log(min=yrange[0],max=yrange[1],
                                           title=yleftlabel, painter=c1y,
                                           parter=None)
                yaxis2 = pyx.graph.axis.log(min=yrange[0],max=yrange[1],
                                            title=yrightlabel, painter=c2y,
                                            parter=None)
        else:
            if yticklabels:
                yaxis = pyx.graph.axis.lin(min=yrange[0],max=yrange[1],
                                           title=yleftlabel, painter=c1y)
                yaxis2 = pyx.graph.axis.lin(min=yrange[0],max=yrange[1],
                                            title=yrightlabel, painter=c2y)
            else:
                yaxis = pyx.graph.axis.lin(min=yrange[0],max=yrange[1],
                                           title=yleftlabel, painter=c1y,
                                           parter=None)
                yaxis2 = pyx.graph.axis.lin(min=yrange[0],max=yrange[1],
                                            title=yrightlabel, painter=c2y,
                                            parter=None)

        if bare_axes:
            if ylog:
                yaxis = pyx.graph.axis.log(min=yrange[0], max=yrange[1],
                                           title=yleftlabel, parter=None)
                yaxis2 = pyx.graph.axis.log(min=yrange[0], max=yrange[1],
                                            title=yrightlabel, parter=None)
            else:
                yaxis = pyx.graph.axis.lin(min=yrange[0], max=yrange[1],
                                           title=yleftlabel, parter=None)
                yaxis2 = pyx.graph.axis.lin(min=yrange[0], max=yrange[1],
                                            title=yrightlabel, parter=None)
            if xlog:
                xaxis = pyx.graph.axis.log(min=xrange[0], max=xrange[1],
                                           title=xleftlabel, parter=None)
                xaxis2 = pyx.graph.axis.log(min=xrange[0], max=xrange[1],
                                            title=xrightlabel, parter=None)
            else:
                xaxis = pyx.graph.axis.lin(min=xrange[0], max=xrange[1],
                                           title=xleftlabel, parter=None)
                xaxis2 = pyx.graph.axis.lin(min=xrange[0], max=xrange[1],
                                            title=xrightlabel, parter=None)

        blank_data = pyx.graph.data.points([(-1e20,-1e20),(-1e19,-1e19)], x=1,y=2)
        if self.canvas is None:
            self.canvas = pyx.graph.graphxy \
                          (width=psize[0], height=psize[1],
                           x=xaxis, y=yaxis, x2=xaxis2, y2=yaxis2,
                           xpos=pos[0], ypos=pos[1])
            if xdata is None:
                self.canvas.plot(blank_data)
            else:
                data = pyx.graph.data.points(np.array([xdata, ydata]).T, x=1, y=2)
                self.canvas.plot(data, [pyx.graph.style.line([pyx.style.linewidth.Thick])])
        else:
            plot = pyx.graph.graphxy \
                   (width=psize[0], height=psize[1],
                    x=xaxis, y=yaxis, x2=xaxis2, y2=yaxis2,
                    xpos=pos[0], ypos=pos[1])
            if xdata is None:
                plot.plot(blank_data)
            else:
                data = pyx.graph.data.points(np.array([xdata, ydata]).T, x=1, y=2)
                plot.plot(data, [pyx.graph.style.line([pyx.style.linewidth.Thick])])
            self.canvas.insert(plot)
        self.axes_drawn = True

#=============================================================================
    
    def axis_box_yt(self, plot, units=None, bare_axes=False,
                    tickcolor=None, xlabel=None, ylabel=None, **kwargs):
        r"""Wrapper around DualEPS.axis_box to automatically fill in the
        axis ranges and labels from a yt plot.

        This also accepts any parameters that DualEPS.axis_box takes.

        Parameters
        ----------
        plot : `yt.visalization.plot_window.PlotWindow`
            yt plot on which the axes are based.
        units : string
            Unit description that overrides yt's unit description.  Only
            affects the axis label.
        bare_axes : boolean
            Set to true to have no annotations or tick marks on all of the
            axes.

        Examples
        --------
        >>> p = SlicePlot(ds, 0, 'density')
        >>> d = DualEPS()
        >>> d.axis_box_yt(p)
        >>> d.save_fig()
        """
       
        if isinstance(plot, (PlotWindow, PhasePlot)):
            plot.refresh()
        if isinstance(plot, PlotWindow):
            data = plot.frb
            width = plot.width[0]
            if units is None:
                units = plot.ds.get_smallest_appropriate_unit(width)
            width = width.in_units(str(units))
            xc = 0.5*(plot.xlim[0] + plot.xlim[1])
            yc = 0.5*(plot.ylim[0] + plot.ylim[1])
            _xrange = [(plot.xlim[i] - xc).in_units(units) for i in (0, 1)]
            _yrange = [(plot.ylim[i] - yc).in_units(units) for i in (0, 1)]
            _xlog = False
            _ylog = False
            if bare_axes:
                _xlabel = ""
                _ylabel = ""
            else:
                if xlabel is not None:
                    _xlabel = xlabel
                else:
                    if data.axis != 4:
                        xi = plot.ds.coordinates.x_axis[data.axis]
                        x_name = plot.ds.coordinates.axis_name[xi]
                        _xlabel = '%s (%s)' % (x_name, units)
                    else:
                        _xlabel = 'x (%s)' % (units)
                if ylabel is not None:
                    _ylabel = ylabel
                else:
                    if data.axis != 4:
                        yi = plot.ds.coordinates.y_axis[data.axis]
                        y_name = plot.ds.coordinates.axis_name[yi]
                        _ylabel = '%s (%s)' % (y_name, units)
                    else:
                        _ylabel = 'y (%s)' % (units)
            if tickcolor is None:
                _tickcolor = pyx.color.cmyk.white
        elif isinstance(plot, ProfilePlot):
            subplot = plot.axes.values()[0]
            # limits for axes
            xlimits = subplot.get_xlim()
            _xrange = (YTQuantity(xlimits[0], 'm'), YTQuantity(xlimits[1], 'm')) # unit hardcoded but afaik it is not used anywhere so it doesn't matter
            if list(plot.axes.ylim.viewvalues())[0][0] is None:
                ylimits = subplot.get_ylim()
            else:
                ylimits = list(plot.axes.ylim.viewvalues())[0]
            _yrange = (YTQuantity(ylimits[0], 'm'), YTQuantity(ylimits[1], 'm')) # unit hardcoded but afaik it is not used anywhere so it doesn't matter
            # axis labels
            xaxis = subplot.xaxis
            _xlabel = pyxize_label(xaxis.label.get_text())
            yaxis = subplot.yaxis
            _ylabel = pyxize_label(yaxis.label.get_text())
            # set log if necessary
            if subplot.get_xscale() == "log":
                 _xlog = True 
            else:
                 _xlog = False
            if subplot.get_yscale() == "log":
                 _ylog = True 
            else:
                 _ylog = False
            _tickcolor = None 
        elif isinstance(plot, PhasePlot):
            k = plot.plots.keys()[0]
            _xrange = plot[k].axes.get_xlim()
            _yrange = plot[k].axes.get_ylim()
            _xlog = plot.profile.x_log
            _ylog = plot.profile.y_log
            if bare_axes:
                _xlabel = ""
                _ylabel = ""
            else:
                if xlabel is not None:
                    _xlabel = xlabel
                else:
                    _xlabel = plot[k].axes.get_xlabel()
                if ylabel is not None:
                    _ylabel = ylabel
                else:
                    _ylabel = plot[k].axes.get_ylabel()
                _xlabel = pyxize_label(_xlabel)
                _ylabel = pyxize_label(_ylabel)
            if tickcolor is None:
                _tickcolor = None
        elif isinstance(plot, np.ndarray):
            ax = plt.gca()
            _xrange = ax.get_xlim()
            _yrange = ax.get_ylim()
            _xlog=False
            _ylog=False
            if bare_axes:
                _xlabel = ""
                _ylabel = ""
            else:
                if xlabel is not None:
                    _xlabel = xlabel
                else:
                    _xlabel = ax.get_xlabel()
                if ylabel is not None:
                    _ylabel = ylabel
                else:
                    _ylabel = ax.get_ylabel()
            if tickcolor is None:
                _tickcolor = None
        else:
            _xrange = plot._axes.get_xlim()
            _yrange = plot._axes.get_ylim()
            _xlog = plot._log_x
            _ylog = plot._log_y
            if bare_axes:
                _xlabel = ""
                _ylabel = ""
            else:
                if xlabel is not None:
                    _xlabel = xlabel
                else:
                    _xlabel = plot._x_label
                if ylabel is not None:
                    _ylabel = ylabel
                else:
                    _ylabel = plot._y_label
            if tickcolor is None:
                _tickcolor = None
        if tickcolor is not None:
            _tickcolor = tickcolor
        self.axis_box(xrange=_xrange, yrange=_yrange, xlabel=_xlabel,
                      ylabel=_ylabel, tickcolor=_tickcolor, xlog=_xlog,
                      ylog=_ylog, bare_axes=bare_axes, **kwargs)

#=============================================================================

    def insert_image(self, filename, pos=(0,0), size=None):
        r"""Inserts a JPEG file in the figure.

        Parameters
        ----------
        filename : string
            Name of the JPEG file
        pos : tuple of floats
            Position of the origin of the image in centimeters
        size : tuple of flots
            Size of image in units of figsize

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.insert_image("image.jpg")
        >>> d.save_fig()
        """
        if size is not None:
            width = size[0]*self.figsize[0]
            height = size[1]*self.figsize[1]
        else:
            width = self.figsize[0]
            height = self.figsize[1]
        image = pyx.bitmap.jpegimage(filename)
        if self.canvas is None:
            self.canvas = pyx.canvas.canvas()
        self.canvas.insert(pyx.bitmap.bitmap(pos[0], pos[1], image,
                                             compressmode=None,
                                             width=width,
                                             height=height))

#=============================================================================

    def insert_image_yt(self, plot, field=None, pos=(0,0), scale=1.0):
        r"""Inserts a bitmap taken from a yt plot.

        Parameters
        ----------
        plot : `yt.visalization.plot_window.PlotWindow`
            yt plot that provides the image
        pos : tuple of floats
            Position of the origin of the image in centimeters.

        Examples
        --------
        >>> p = SlicePlot(ds, 0, 'density')
        >>> d = DualEPS()
        >>> d.axis_box_yt(p)
        >>> d.insert_image_yt(p)
        >>> d.save_fig()

        Notes
        -----
        For best results, set use_colorbar=False when creating the yt
        image.
        """
        from ._mpl_imports import FigureCanvasAgg

        # We need to remove the colorbar (if necessary), remove the
        # axes, and resize the figure to span the entire figure
        force_square = False
        if self.canvas is None:
            self.canvas = pyx.canvas.canvas()
        if isinstance(plot, (PlotWindow, PhasePlot)):
            if field is None:
                self.field = plot.plots.keys()[0]
                mylog.warning("No field specified.  Choosing first field (%s)" % \
                              str(self.field))
            else:
                self.field = plot.data_source._determine_fields(field)[0]
            if self.field not in plot.plots.keys():
                raise RuntimeError("Field '%s' does not exist!" % str(self.field))
            if isinstance(plot, PlotWindow):
                plot.hide_colorbar()
                plot.hide_axes()
            else:
                plot.plots[self.field]._toggle_axes(False)
                plot.plots[self.field]._toggle_colorbar(False)
            plot.refresh()
            _p1 = plot.plots[self.field].figure
            force_square = True
        elif isinstance(plot, ProfilePlot):
            _p1 = plot.figures.items()[0][1]
        elif isinstance(plot, np.ndarray):
            plt.figure()
            iplot = plt.figimage(plot)
            _p1 = iplot.figure
            _p1.set_size_inches(self.figsize[0], self.figsize[1])
            ax = plt.gca()
            _p1.add_axes(ax)
        else:
            raise RuntimeError("Unknown plot type")

        _p1.axes[0].set_position([0,0,1,1])  # rescale figure
        _p1.set_facecolor('w')  # set background color
        figure_canvas = FigureCanvasAgg(_p1)
        figure_canvas.draw()
        size = (_p1.get_size_inches() * _p1.dpi).astype('int')

        # Account for non-square images after removing the colorbar.
        scale *= 1.0 - 1.0 / (_p1.dpi * self.figsize[0])
        if force_square:
            yscale = scale * float(size[1]) / float(size[0])
        else:
            yscale = scale
        image = pyx.bitmap.image(size[0], size[1], "RGB",
                                 figure_canvas.tostring_rgb())
        self.canvas.insert(pyx.bitmap.bitmap(pos[0], pos[1], image,
                                             width=scale*self.figsize[0],
                                             height=yscale*self.figsize[1]))

#=============================================================================

    def colorbar(self, name, zrange=(0,1), label="", log=False, tickcolor=None,
                 orientation="right", pos=[0,0], shrink=1.0):
        r"""Places a colorbar adjacent to the current figure.

        Parameters
        ----------
        name : string
            name of the (matplotlib) colormap to use
        zrange : tuple of floats
            min and max of the colorbar's range
        label : string
            colorbar label
        log : boolean
            Flag to use a logarithmic scale
        tickcolor : `pyx.color.*.*`
            Color for the tickmarks.  Example: pyx.color.cmyk.black
        orientation : string
            Placement of the colorbar.  Can be "left", "right", "top",
            or "bottom".
        pos : list of floats
            (x,y) position of the origin of the colorbar in centimeters.
        shrink : float
            Factor to shrink the colorbar's size.  A value of 1 means the
            colorbar will have a height / width of the figure.

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.insert_image("image.jpg")
        >>> d.colorbar("hot", xrange=(1e-2, 1e-4), log=True,
                       label="Density [cm$^{-3}$]")
        >>> d.save_fig()
        """
        if orientation == "right":
            origin = (pos[0]+self.figsize[0]+0.5, pos[1])
            size = (0.1*self.figsize[0], self.figsize[1])
            imsize = (1,256)
        elif orientation == "left":
            origin = (pos[0]-0.5-0.1*self.figsize[0], pos[1])
            size = (0.1*self.figsize[0], self.figsize[1])
            imsize = (1,256)
        elif orientation == "top":
            origin = (pos[0], pos[1]+self.figsize[1]+0.5)
            imorigin = (pos[0]+self.figsize[0], pos[1]+self.figsize[1]+0.5)
            size = (self.figsize[0], 0.1*self.figsize[1])
            imsize = (256,1)
        elif orientation == "bottom":
            origin = (pos[0], pos[1]-0.5-0.1*self.figsize[1])
            imorigin = (pos[0]+self.figsize[0], pos[1]-0.5-0.1*self.figsize[1])
            size = (self.figsize[0], 0.1*self.figsize[1])
            imsize = (256,1)
        else:
            raise RuntimeError("orientation %s unknown" % orientation)
            return

        # If shrink is a scalar, then convert into tuple
        if not isinstance(shrink, (tuple,list)):
            shrink = (shrink, shrink)

        # Scale the colorbar
        shift = (0.5*(1.0-shrink[0])*size[0], 0.5*(1.0-shrink[1])*size[1])
        # To facilitate strething rather than shrinking
        # If stretched in both directions (makes no sense?) then y dominates. 
        if(shrink[0] > 1.0):
            shift = (0.05*self.figsize[0], 0.5*(1.0-shrink[1])*size[1])
        if(shrink[1] > 1.0):
            shift = (0.5*(1.0-shrink[0])*size[0], 0.05*self.figsize[1])
        size = (size[0] * shrink[0], size[1] * shrink[1])
        origin = (origin[0] + shift[0], origin[1] + shift[1])

        # Convert the colormap into a string
        x = np.linspace(1,0,256)
        cm_string = cm.cmap_d[name](x, bytes=True)[:,0:3].tostring()

        cmap_im = pyx.bitmap.image(imsize[0], imsize[1], "RGB", cm_string)
        if orientation == "top" or orientation == "bottom":
            imorigin = (imorigin[0] - shift[0], imorigin[1] + shift[1])
            self.canvas.insert(pyx.bitmap.bitmap(imorigin[0], imorigin[1], cmap_im,
                                                 width=-size[0], height=size[1]))
        else:
            self.canvas.insert(pyx.bitmap.bitmap(origin[0], origin[1], cmap_im,
                                                 width=size[0], height=size[1]))

        if tickcolor is None:
            c1 = pyx.graph.axis.painter.regular(tickattrs=[pyx.color.cmyk.black])
            pyx.graph.axis.painter.regular(tickattrs=[pyx.color.cmyk.black],
                                           labelattrs=None)
        else:
            c1 = pyx.graph.axis.painter.regular(tickattrs=[tickcolor])
            pyx.graph.axis.painter.regular(tickattrs=[tickcolor],
                                           labelattrs=None)
        if log:
            yaxis = pyx.graph.axis.log(min=zrange[0],max=zrange[1],
                                       title=label, painter=c1)
            yaxis2 = pyx.graph.axis.log(min=zrange[0],max=zrange[1],parter=None)
        else:
            yaxis = pyx.graph.axis.lin(min=zrange[0],max=zrange[1],
                                       title=label, painter=c1)
            yaxis2 = pyx.graph.axis.lin(min=zrange[0], max=zrange[1], parter=None)
        xaxis = pyx.graph.axis.lin(parter=None)

        if orientation == "right":
            _colorbar = pyx.graph.graphxy(width=size[0], height=size[1],
                                          xpos=origin[0], ypos=origin[1],
                                          x=xaxis, y=yaxis2, y2=yaxis)
        elif orientation == "left":
            _colorbar = pyx.graph.graphxy(width=size[0], height=size[1],
                                          xpos=origin[0], ypos=origin[1],
                                          x=xaxis, y2=yaxis2, y=yaxis)
        elif orientation == "top":
            _colorbar = pyx.graph.graphxy(width=size[0], height=size[1],
                                          xpos=origin[0], ypos=origin[1],
                                          y=xaxis, x=yaxis2, x2=yaxis)
        elif orientation == "bottom":
            _colorbar = pyx.graph.graphxy(width=size[0], height=size[1],
                                          xpos=origin[0], ypos=origin[1],
                                          y=xaxis, x2=yaxis2, x=yaxis)
            
        
        blank_data = pyx.graph.data.points([(-1e10,-1e10),(-9e10,-9e10)],
                                           x=1, y=2)
        _colorbar.plot(blank_data)
        self.canvas.insert(_colorbar)        

#=============================================================================

    def colorbar_yt(self, plot, field=None, cb_labels = None, **kwargs):
        r"""Wrapper around DualEPS.colorbar to take information from a yt plot.

        Accepts all parameters that DualEPS.colorbar takes.

        Parameters
        ----------
        plot : A yt plot
            yt plot from which the information is taken.
        cb_labels : list of labels for the colorbars. List should be the same
                    size as the number of colorbars used. Should be passed 
                    into this function by either the singleplot or multiplot api.

        Examples
        --------
        >>> p = SlicePlot(ds, 0, 'density')
        >>> p.hide_colorbar()
        >>> d = DualEPS()
        >>> d.axis_box_yt(p)
        >>> d.insert_image_yt(p)
        >>> d.colorbar_yt(p)
        >>> d.save_fig()
        """
        
        if isinstance(plot, ProfilePlot):
            raise RuntimeError("When using ProfilePlots you must either set yt_nocbar=True or provide colorbar flags so that the profiles don't have colorbars")
        _cmap = None
        if field is not None:
            self.field = plot.data_source._determine_fields(field)[0]
        if isinstance(plot, (PlotWindow, PhasePlot)):
            _cmap = plot._colormaps[self.field]
        else:
            if plot.cmap is not None:
                _cmap = plot.cmap.name
        if _cmap is None:
            _cmap = ytcfg.get("yt", "default_colormap")
        if isinstance(plot, (PlotWindow, PhasePlot)):
            if isinstance(plot, PlotWindow):
                try:
                    _zlabel = plot.frb[self.field].info["label"]
                    _unit = Unit(plot.frb[self.field].units, 
                                 registry=plot.ds.unit_registry)
                    units = _unit.latex_representation()
                    # PyX does not support \frac because it's based on TeX.
                    units = pyxize_label(units)
                    _zlabel += r' (' + units + r')'
                except NotImplementedError: 
                    print("Colorbar label not available")
                    _zlabel = ''
            else:
                _, _, z_title = plot._get_field_title(self.field, plot.profile)
                _zlabel = pyxize_label(z_title)
            _zlabel = _zlabel.replace("_","\;")
            _zlog = plot.get_log(self.field)[self.field]
            if plot.plots[self.field].zmin is None:
                zmin = plot.plots[self.field].image._A.min()
            else:
                zmin = plot.plots[self.field].zmin
            if plot.plots[self.field].zmax is None:
                zmax = plot.plots[self.field].image._A.max()
            else:
                zmax = plot.plots[self.field].zmax
            _zrange = (zmin, zmax)
        else:
            _zlabel = plot._z_label.replace("_","\;")
            _zlog = plot._log_z
            _zrange = (plot.norm.vmin, plot.norm.vmax)
        if cb_labels is not None:  #Overrides deduced labels
            _zlabel = cb_labels.pop()
        self.colorbar(_cmap, zrange=_zrange, label=_zlabel, log=_zlog, **kwargs)

#=============================================================================

    def circle(self, radius=0.2, loc=(0.5,0.5),
               color=pyx.color.cmyk.white,
               linewidth=pyx.style.linewidth.normal):
        r"""Draws a circle in the current figure.

        Parameters
        ----------
        radius : float
            Radius of the circle in units of figsize
        loc : tuple of floats
            Location of the circle's center in units of figsize
        color : `pyx.color.*.*`
            Color of the circle stroke.  Example: pyx.color.cmyk.white
        linewidth : `pyx.style.linewidth.*`
            Width of the circle stroke width. Example:
            pyx.style.linewidth.normal

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.insert_image("image.jpg")
        >>> d.circle(radius=0.1, color=pyx.color.cmyk.Red)
        >>> d.save_fig()
        """
        circle = pyx.path.circle(self.figsize[0]*loc[0],
                                 self.figsize[1]*loc[1],
                                 self.figsize[0]*radius)
        self.canvas.stroke(circle, [color, linewidth])

#=============================================================================

    def arrow(self, size=0.2, label="", loc=(0.05,0.08), labelloc="top",
              color=pyx.color.cmyk.white,
              linewidth=pyx.style.linewidth.normal):
        r"""Draws an arrow in the current figure

        Parameters
        ----------
        size : float
            Length of arrow (base to tip) in units of the figure size.
        label : string
            Annotation label of the arrow.
        loc : tuple of floats
            Location of the left hand side of the arrow in units of
            the figure size.
        labelloc : string
            Location of the label with respect to the line.  Can be
            "top" or "bottom"
        color : `pyx.color.*.*`
            Color of the arrow.  Example: pyx.color.cymk.white
        linewidth : `pyx.style.linewidth.*`
            Width of the arrow.  Example: pyx.style.linewidth.normal

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.insert_image("arrow_image.jpg")
        >>> d.arrow(size=0.2, label="Black Hole!", loc=(0.05, 0.1))
        >>> d.save_fig()
        """
        line = pyx.path.line(self.figsize[0]*loc[0],
                             self.figsize[1]*loc[1],
                             self.figsize[0]*(loc[0]+size),
                             self.figsize[1]*loc[1])
        self.canvas.stroke(line, [linewidth, color, pyx.deco.earrow()])
       

        if labelloc == "bottom":
            yoff = -0.1*size
            valign = pyx.text.valign.top
        else:
            yoff = +0.1*size
            valign = pyx.text.valign.bottom
        if label != "":
            self.canvas.text(self.figsize[0]*(loc[0]+0.5*size),
                             self.figsize[1]*(loc[1]+yoff), label,
                             [color, valign, pyx.text.halign.center])

        


#=============================================================================

    def scale_line(self, size=0.2, label="", loc=(0.05,0.08), labelloc="top",
                   color=pyx.color.cmyk.white,
                   linewidth=pyx.style.linewidth.normal):
        r"""Draws a scale line in the current figure.

        Parameters
        ----------
        size : float
            Length of the scale line in units of the figure size.
        label : string
            Annotation label of the scale line.
        loc : tuple of floats
            Location of the left hand side of the scale line in units of
            the figure size.
        labelloc : string
            Location of the label with respect to the line.  Can be
            "top" or "bottom"
        color : `pyx.color.*.*`
            Color of the scale line.  Example: pyx.color.cymk.white
        linewidth : `pyx.style.linewidth.*`
            Width of the scale line.  Example: pyx.style.linewidth.normal

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.insert_image("image.jpg")
        >>> d.scale_line(size=0.2, label="1 kpc", loc=(0.05, 0.1))
        >>> d.save_fig()
        """
        
        line = pyx.path.line(self.figsize[0]*loc[0],
                             self.figsize[1]*loc[1],
                             self.figsize[0]*(loc[0]+size),
                             self.figsize[1]*loc[1])
        self.canvas.stroke(line, [linewidth, color])
        line = pyx.path.line(self.figsize[0]*loc[0],
                             self.figsize[1]*(loc[1]-0.1*size),
                             self.figsize[0]*loc[0],
                             self.figsize[1]*(loc[1]+0.1*size))
        self.canvas.stroke(line, [linewidth, color])
        line = pyx.path.line(self.figsize[0]*(loc[0]+size),
                             self.figsize[1]*(loc[1]-0.1*size),
                             self.figsize[0]*(loc[0]+size),
                             self.figsize[1]*(loc[1]+0.1*size))
        self.canvas.stroke(line, [linewidth, color])

        if labelloc == "bottom":
            yoff = -0.1*size
            valign = pyx.text.valign.top
        else:
            yoff = +0.1*size
            valign = pyx.text.valign.bottom
        if label != "":
            self.canvas.text(self.figsize[0]*(loc[0]+0.5*size),
                             self.figsize[1]*(loc[1]+yoff), label,
                             [color, valign, pyx.text.halign.center])

#=============================================================================

    def title_box(self, text, color=pyx.color.cmyk.black,
                  bgcolor=pyx.color.cmyk.white, loc=(0.02,0.98),
                  halign=pyx.text.halign.left,
                  valign=pyx.text.valign.top,
                  text_opts=[]):
        r"""Inserts a box with text in the current figure.

        Parameters
        ----------
        text : string
            String to insert in the textbox.
        color : `pyx.color.*.*`
            Color of the text.  Example: pyx.color.cmyk.black
        bgcolor : `pyx.color.*.*`
            Color of the textbox background.  Example: pyx.color.cmyk.white
        loc : tuple of floats
            Location of the textbox origin in units of the figure size.
        halign : `pyx.text.halign.*`
            Horizontal alignment of the text.  Example: pyx.text.halign.left
        valign : `pyx.text.valign.*`
            Vertical alignment of the text.  Example: pyx.text.valign.top

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.insert_image("image.jpg")
        >>> d.title_box("Halo 1", loc=(0.05,0.95))
        >>> d.save_fig()
        """
        tbox = self.canvas.text(self.figsize[0]*loc[0],
                                self.figsize[1]*loc[1],
                                text, [color, valign, halign] + text_opts)
        if bgcolor is not None:
            tpath = tbox.bbox().enlarged(2*pyx.unit.x_pt).path()
            self.canvas.draw(tpath, [pyx.deco.filled([bgcolor]),
                                     pyx.deco.stroked()])
        self.canvas.insert(tbox)
        
#=============================================================================

    def save_fig(self, filename="test", format="eps", resolution=250):
        r"""Saves current figure to a file.

        Parameters
        ----------
        filename : string
            Name of the saved file without the extension.
        format : string
            Format type.  Can be "eps" or "pdf"

        Examples
        --------
        >>> d = DualEPS()
        >>> d.axis_box(xrange=(0,100), yrange=(1e-3,1), ylog=True)
        >>> d.save_fig("image1", format="pdf")
        """
        if format =="eps":
            self.canvas.writeEPSfile(filename)
        elif format == "pdf":
            self.canvas.writePDFfile(filename)
        elif format == "png":
             self.canvas.writeGSfile(filename+".png", "png16m", resolution=resolution)
        elif format == "jpg":
             self.canvas.writeGSfile(filename+".jpeg", "jpeg", resolution=resolution)
        else:
            raise RuntimeError("format %s unknown." % (format))
            
#=============================================================================
#=============================================================================
#=============================================================================

def multiplot(ncol, nrow, yt_plots=None, fields=None, images=None, 
              xranges=None, yranges=None, xlabels=None, ylabels=None,
              xdata=None, ydata=None, colorbars=None,
              shrink_cb=0.95, figsize=(8,8), margins=(0,0), titles=None,
              savefig=None, format="eps", yt_nocbar=False, bare_axes=False,
              xaxis_flags=None, yaxis_flags=None,
              cb_flags=None, cb_location=None, cb_labels=None):
    r"""Convenience routine to create a multi-panel figure from yt plots or
    JPEGs.  The images are first placed from the origin, and then
    bottom-to-top and left-to-right.

    Parameters
    ----------
    ncol : integer
        Number of columns in the figure.
    nrow : integer
        Number of rows in the figure.
    yt_plots : list of yt plot instances
        yt plots to include in the figure.
    images : list of strings
        JPEG filenames to include in the figure.
    xranges : list of tuples
        The min and max of the x-axes
    yranges : list of tuples
        The min and max of the y-axes
    xlabels : list of strings
        Labels for the x-axes
    ylabels : list of strings
        Labels for the y-axes
    colorbars : list of dicts
        Dicts that describe the type of colorbar to be used in each
        figure.  Use the function return_cmap() to create these dicts.
    shrink_cb : float
        Factor by which the colorbar is shrunk.
    figsize : tuple of floats
        The width and height of a single figure in centimeters.
    margins : tuple of floats
        The horizontal and vertical margins between panels in centimeters.
    titles : list of strings
        Titles that are placed in textboxes in each panel.
    savefig : string
        Name of the saved file without the extension.
    format : string
        File format of the figure. eps or pdf accepted.
    yt_nocbar : boolean
        Flag to indicate whether or not colorbars are created.
    bare_axes : boolean
        Set to true to have no annotations or tick marks on all of the
        axes.
    cb_flags : list of booleans
        Flags for each plot to have a colorbar or not.
    cb_location : list of strings
        Strings to control the location of the colorbar (left, right, 
        top, bottom)
    cb_labels : list of labels for the colorbars. List should be the same
                size as the number of colorbars used.

    Examples
    --------
    >>> images = ["density.jpg", "hi_density.jpg", "entropy.jpg",
    >>>           "special.jpg"]
    >>> cbs=[]
    >>> cbs.append(return_cmap("arbre", "Density [cm$^{-3}$]", (0,10), False))
    >>> cbs.append(return_cmap("kelp", "HI Density", (0,5), False))
    >>> cbs.append(return_cmap("hot", r"Entropy [K cm$^2$]", (1e-2,1e6), True))
    >>> cbs.append(return_cmap("Spectral", "Stuff$_x$!", (1,300), True))
    >>> 
    >>> mp = multiplot(2,2, images=images, margins=(0.1,0.1),
    >>>                titles=["1","2","3","4"],
    >>>                xlabels=["one","two"], ylabels=None, colorbars=cbs,
    >>>                shrink_cb=0.95)
    >>> mp.scale_line(label="$r_{vir}$", labelloc="top")
    >>> mp.save_fig("multiplot")

    Notes
    -----
    If given both yt_plots and images, this will get preference to the
    yt plots.
    """
    # Error check
    npanels = ncol*nrow
    if(cb_labels is not None):
        cb_labels.reverse()   #Because I pop the list
    
    if images is not None:
        if len(images) != npanels:
            raise RuntimeError("Number of images (%d) doesn't match nrow(%d)"\
                               " x ncol(%d)." % (len(images), nrow, ncol))
            return
    if yt_plots is None and images is None:
        raise RuntimeError("Must supply either yt_plots or image filenames.")
        return
    if yt_plots is not None and images is not None:
        mylog.warning("Given both images and yt plots.  Ignoring images.")
    if yt_plots is not None:
        _yt = True
    else:
        _yt = False
    if fields is None:
        fields = [None] * npanels

    # If no ranges or labels given and given only images, fill them in.
    if not _yt:
        if xranges is None:
            xranges = []
            for i in range(npanels): xranges.append((0,1))
        if yranges is None:
            yranges = []
            for i in range(npanels): yranges.append((0,1))
        if xlabels is None:
            xlabels = []
            for i in range(npanels): xlabels.append("")
        if ylabels is None:
            ylabels = []
            for i in range(npanels): ylabels.append("")

    d = DualEPS(figsize=figsize)
    for j in range(nrow):
        invj = nrow - j - 1
        ypos = invj*(figsize[1] + margins[1])
        for i in range(ncol):
            xpos = i*(figsize[0] + margins[0])
            index = j*ncol + i
            if isinstance(yt_plots, list):
                this_plot = yt_plots[index]
            else:
                this_plot = yt_plots
            if j == nrow-1:
                xaxis = 0
            elif j == 0:
                xaxis = 1
            else:
                xaxis = -1
            if i == 0:
                yaxis = 0
            elif i == ncol-1:
                yaxis = 1
            else:
                yaxis = -1
            if xdata is None:
                _xdata = None
            else:
                _xdata = xdata[index]
            if ydata is None:
                _ydata = None
            else:
                _ydata = ydata[index]
            if xaxis_flags is not None:
                if xaxis_flags[index] is not None:
                    xaxis = xaxis_flags[index]
            if yaxis_flags is not None:
                if yaxis_flags[index] is not None:
                    yaxis = yaxis_flags[index]
            if _yt:
                this_plot._setup_plots()
                if xlabels is not None:
                    xlabel = xlabels[i]
                else:
                    xlabel = None
                if ylabels is not None:
                    ylabel = ylabels[j]
                else:
                    ylabel = None
                d.insert_image_yt(this_plot, pos=(xpos, ypos),
                                  field=fields[index])
                d.axis_box_yt(this_plot, pos=(xpos, ypos),
                              bare_axes=bare_axes, xaxis_side=xaxis,
                              yaxis_side=yaxis,
                              xlabel=xlabel, ylabel=ylabel,
                              xdata=_xdata, ydata=_ydata)
            else:
                d.insert_image(images[index], pos=(xpos,ypos))
                d.axis_box(pos = (xpos, ypos),
                           xrange=xranges[index], yrange=yranges[index],
                           xlabel=xlabels[i], ylabel=ylabels[j],
                           bare_axes=bare_axes, xaxis_side=xaxis, yaxis_side=yaxis,
                           xdata=_xdata, ydata=_ydata)
            if titles is not None:
                if titles[index] is not None:
                    d.title_box(titles[index],
                                loc=(i+0.05+i*margins[0]/figsize[0],
                                     j+0.98+j*margins[1]/figsize[1]))

    # Insert colorbars after all axes are placed because we want to
    # put them on the edges of the bounding box.
    bbox = (100.0 * d.canvas.bbox().left().t,
            100.0 * d.canvas.bbox().right().t - d.figsize[0],
            100.0 * d.canvas.bbox().bottom().t,
            100.0 * d.canvas.bbox().top().t - d.figsize[1])
    for j in range(nrow):
        invj = nrow - j - 1
        ypos0 = invj*(figsize[1] + margins[1])
        for i in range(ncol):
            xpos0 = i*(figsize[0] + margins[0])
            index = j*ncol + i
            if isinstance(yt_plots, list):
                this_plot = yt_plots[index]
            else:
                this_plot = yt_plots
            if (not _yt and colorbars is not None) or (_yt and not yt_nocbar):
                if cb_flags is not None:
                    if not cb_flags[index]:
                        continue
                if cb_location is None:
                    if ncol == 1:
                        orientation = "right"
                    elif i == 0:
                        orientation = "left"
                    elif i+1 == ncol:
                        orientation = "right"
                    elif j == 0:
                        orientation = "bottom"
                    elif j+1 == nrow:
                        orientation = "top"
                    else:
                        orientation = None  # Marker for interior plot
                else:
                    if isinstance(cb_location, dict):
                        if fields[index] not in cb_location.keys():
                            raise RuntimeError("%s not found in cb_location dict" %
                                               fields[index])
                            return
                        orientation = cb_location[fields[index]]
                    elif isinstance(cb_location, list):
                        orientation = cb_location[index]
                    else:
                        raise RuntimeError("Bad format: cb_location")
                if orientation == "right":
                    xpos = bbox[1]
                    ypos = ypos0
                elif orientation == "left":
                    xpos = bbox[0]
                    ypos = ypos0
                elif orientation == "bottom":
                    ypos = bbox[2]
                    xpos = xpos0
                elif orientation == "top":
                    ypos = bbox[3]
                    xpos = xpos0
                else:
                    mylog.warning("Unknown colorbar location %s. "
                                  "No colorbar displayed." % orientation)
                    orientation = None  # Marker for interior plot

                if orientation is not None:
                    if _yt:
                        # Set field if undefined
                        if fields[index] is None:
                            fields[index] = d.return_field(yt_plots[index])
                                              
                        d.colorbar_yt(this_plot,
                                      field=fields[index],
                                      pos=[xpos,ypos],
                                      shrink=shrink_cb,
                                      orientation=orientation,
                                      cb_labels=cb_labels)
                    else:
                        d.colorbar(colorbars[index]["cmap"],
                                   zrange=colorbars[index]["range"],
                                   label=colorbars[index]["name"],
                                   log=colorbars[index]["log"],
                                   orientation=orientation,
                                   pos=[xpos,ypos],
                                   shrink=shrink_cb)

    if savefig is not None:
        d.save_fig(savefig, format=format)

    return d

#=============================================================================

def multiplot_yt(ncol, nrow, plots, fields=None, **kwargs):
    r"""Wrapper for multiplot that takes a yt PlotWindow

    Accepts all parameters used in multiplot.

    Parameters
    ----------
    ncol : integer
        Number of columns in the figure.
    nrow : integer
        Number of rows in the figure.
    plots : ``PlotWindow`` instance, ``PhasePlot`` instance, or list of plots
        yt plots to be used.

    Examples
    --------
    >>> p1 = SlicePlot(ds, 0, 'density')
    >>> p1.set_width(10, 'kpc')
    >>>
    >>> p2 = SlicePlot(ds, 0, 'temperature')
    >>> p2.set_width(10, 'kpc')
    >>> p2.set_cmap('temperature', 'hot')
    >>>
    >>> sph = ds.sphere(ds.domain_center, (10, 'kpc'))
    >>> p3 = PhasePlot(sph, 'radius', 'density', 'temperature',
    ...                weight_field='cell_mass')
    >>>
    >>> p4 = PhasePlot(sph, 'radius', 'density', 'pressure', 'cell_mass')
    >>>
    >>> mp = multiplot_yt(2, 2, [p1, p2, p3, p4], savefig="yt", shrink_cb=0.9,
    ...                   bare_axes=True, yt_nocbar=False, margins=(0.5,0.5))
    """
    # Determine whether the plots are organized in a PlotWindow, or list
    # of PlotWindows
    if isinstance(plots, (PlotWindow, PhasePlot)):
        if fields is None:
            fields = plots.fields
        if len(fields) < nrow*ncol:
            raise RuntimeError("Number of plots is less "\
                               "than nrow(%d) x ncol(%d)." % \
                               (len(fields), nrow, ncol))
            return
        figure = multiplot(ncol, nrow, yt_plots=plots, fields=fields, **kwargs)
    elif isinstance(plots, list) and isinstance(plots[0], (PlotWindow, PhasePlot)):
        if len(plots) < nrow*ncol:
            raise RuntimeError("Number of plots is less "\
                               "than nrow(%d) x ncol(%d)." % \
                               (len(fields), nrow, ncol))
            return
        figure = multiplot(ncol, nrow, yt_plots=plots, fields=fields, **kwargs)
    else:
        raise RuntimeError("Unknown plot type in multiplot_yt")
        return
    return figure

#=============================================================================

def single_plot(plot, field=None, figsize=(12,12), cb_orient="right", 
                bare_axes=False, savefig=None, colorbar=True, 
                file_format='eps', **kwargs):
    r"""Wrapper for DualEPS routines to create a figure directy from a yt
    plot.  Calls insert_image_yt, axis_box_yt, and colorbar_yt.

    Parameters
    ----------
    plot : `yt.visalization.plot_window.PlotWindow`
        yt plot that provides the image and metadata
    figsize : tuple of floats
        Size of the figure in centimeters.
    cb_orient : string
        Placement of the colorbar.  Can be "left", "right", "top", or
        "bottom".
    bare_axes : boolean
        Set to true to have no annotations or tick marks on all of the axes.
    savefig : string
        Name of the saved file without the extension.
    colorbar : boolean
        Set to true to include a colorbar
    file_format : string
        Format type.  Can be "eps" or "pdf"

    Examples
    --------
    >>> p = SlicePlot(ds, 0, 'density')
    >>> p.set_width(0.1,'kpc')
    >>> single_plot(p, savefig="figure1")
    """
    d = DualEPS(figsize=figsize)
    d.insert_image_yt(plot, field=field)
    d.axis_box_yt(plot, bare_axes=bare_axes, **kwargs)
    if colorbar and not isinstance(plot, ProfilePlot):
        d.colorbar_yt(plot, orientation=cb_orient)
    if savefig is not None:
        d.save_fig(savefig, format=file_format)
    return d

#=============================================================================
def return_cmap(cmap=None, label="", range=(0,1), log=False):
    r"""Returns a dict that describes a colorbar.  Exclusively for use with
    multiplot.

    Parameters
    ----------
    cmap : string
        name of the (matplotlib) colormap to use
    label : string
        colorbar label
    range : tuple of floats
        min and max of the colorbar's range
    log : boolean
        Flag to use a logarithmic scale

    Examples
    --------
    >>> cb = return_cmap("arbre", "Density [cm$^{-3}$]", (0,10), False)
    """
    if cmap is None:
        cmap = ytcfg.get("yt", "default_colormap")
    return {'cmap': cmap, 'name': label, 'range': range, 'log': log}
    

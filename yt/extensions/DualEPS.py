"""
DualEPS: A class to combine bitmap compression and vector graphics

Author: John Wise <jwise@astro.princeton.edu>
Date: April 2010
Affiliation: Princeton
Homepage: http://yt.enzotools.org/

Requirements: PyX

License:
  Copyright (C) 2010 John Wise.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import pyx
from yt.mods import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib import cm

class DualEPS:
    def __init__(self, figsize=(12,12)):
        """
        Initializes the DualEPS class with a figure size of *figsize*
        centimeters for a single plot.
        """
        pyx.unit.set(xscale=1.4)
        self.figsize = figsize
        self.canvas = None
        self.colormaps = None
        self.axes_drawn = False

    def hello_world(self):
        """
        A simple test.
        """
        if self.canvas is None:
            self.canvas = pyx.canvas.canvas()
        p = pyx.path.line(0,0,1,1)
        self.canvas.stroke(p)
        self.canvas.text(0,0, "Hello world.")

#=============================================================================

    def axis_box(self, xrange=(0,1), yrange=(0,1), xlabel="", ylabel="",
                 xlog=False, ylog=False, tickcolor=None, bare_axes=False,
                 pos=(0,0), xaxis_side=0, yaxis_side=0):
        """
        Draws an axis box at position *pos* with a range of *xrange*
        and *yrange*, labels *xlabel* and *ylabel*.  The colors of the
        tickmarks can be specified with *tickcolor*.  For no tick
        labels or marks, set *bare_axes* to True.  For the x-axis
        (y-axis) labels to be on the right (top), set *xaxis_side*
        (*yaxis_side*) to 1.
        """
        if tickcolor is None:
            c1 = pyx.graph.axis.painter.regular\
                 (tickattrs=[pyx.color.cmyk.black])
            c2 = pyx.graph.axis.painter.regular\
                 (tickattrs=[pyx.color.cmyk.black], labelattrs=None)
        else:
            c1 = pyx.graph.axis.painter.regular(tickattrs=[tickcolor])
            c2 = pyx.graph.axis.painter.regular\
                 (tickattrs=[tickcolor], labelattrs=None)

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

        blank_data = pyx.graph.data.points([(-1,-1),(-0.99,-0.99)], x=1,y=2)
        if self.canvas is None:
            self.canvas = pyx.graph.graphxy \
                          (width=self.figsize[0], height=self.figsize[1],
                           x=xaxis, y=yaxis, x2=xaxis2, y2=yaxis2,
                           xpos=pos[0], ypos=pos[1])
            self.canvas.plot(blank_data)
        else:
            plot = pyx.graph.graphxy \
                   (width=self.figsize[0], height=self.figsize[1],
                    x=xaxis, y=yaxis, x2=xaxis2, y2=yaxis2,
                    xpos=pos[0], ypos=pos[1])
            plot.plot(blank_data)
            self.canvas.insert(plot)
        self.axes_drawn = True

#=============================================================================

    def axis_box_yt(self, plot, units=None, bare_axes=False, **kwargs):
        plot._redraw_image()
        if isinstance(plot, raven.PlotTypes.VMPlot):
            if units == None:
                # Determine the best units
                astro_units = ['cm', 'rsun', 'au', 'pc', 'kpc', 'Mpc']
                best_fit = 0
                while plot.width*plot.pf[astro_units[best_fit]] > 1e3 and \
                          best_fit < len(astro_units):
                    best_fit += 1
                units = astro_units[best_fit]
            _xrange = (0, plot.width * plot.pf[units])
            _yrange = (0, plot.width * plot.pf[units])
            _xlog = False
            _ylog = False
            if bare_axes:
                _xlabel = ""
                _ylabel = ""
            else:
                _xlabel = '%s (%s)' % (lagos.x_names[plot.data.axis], units)
                _ylabel = '%s (%s)' % (lagos.y_names[plot.data.axis], units)
            _tickcolor = pyx.color.cmyk.white
        else:
            _xrange = plot._axes.get_xlim()
            _yrange = plot._axes.get_ylim()
            _xlog = plot._log_x
            _ylog = plot._log_y
            if bare_axes:
                _xlabel = ""
                _ylabel = ""
            else:
                _xlabel = plot._x_label
                _ylabel = plot._y_label
            _tickcolor = None
        self.axis_box(xrange=_xrange, yrange=_yrange, xlabel=_xlabel,
                      ylabel=_ylabel, tickcolor=_tickcolor, xlog=_xlog,
                      ylog=_ylog, bare_axes=bare_axes, **kwargs)

#=============================================================================

    def insert_image(self, filename, pos=(0,0)):
        """
        Inserts a JPEG file *filename* at *pos*.
        """
        image = pyx.bitmap.jpegimage(filename)
        self.canvas.insert(pyx.bitmap.bitmap(pos[0], pos[1], image,
                                             compressmode=None,
                                             width=self.figsize[0],
                                             height=self.figsize[1]))

#=============================================================================

    def insert_image_yt(self, plot, pos=(0,0)):
        """
        Inserts a bitmap taken from a yt plot.  For best results, set
        use_colorbar=False when creating the yt image.
        """
        # We need to remove the colorbar (if necessary), remove the
        # axes, and resize the figure to span the entire figure
        if plot.colorbar != None and \
               isinstance(plot, raven.PlotTypes.VMPlot):
            print "WARNING: Image (slices, projections, etc.) plots must not"\
                  "have a colorbar."
            print "Removing it."
            plot.colorbar = None
        if self.canvas is None:
            self.canvas = pyx.canvas.canvas()
        plot._redraw_image()
        _p1 = plot._figure
        if isinstance(plot, raven.PlotTypes.ProfilePlot):
            # Remove colorbar
            _p1.delaxes(_p1.axes[1])
        _p1.axes[0].set_axis_off()  # remove axes
        _p1.axes[0].set_position([0,0,1,1])  # rescale figure
        _p1.set_facecolor('w')  # set background color
        figure_canvas = FigureCanvas(_p1)
        figure_canvas.draw()
        size = _p1.get_size_inches() * _p1.dpi
        image = pyx.bitmap.image(size[0], size[1], "RGB",
                                 figure_canvas.tostring_rgb())
        figure_canvas.print_png('test.png')
        self.canvas.insert(pyx.bitmap.bitmap(pos[0], pos[1], image,
                                             width=self.figsize[0],
                                             height=self.figsize[1]))

#=============================================================================

    def colorbar(self, name, zrange=(0,1), label="", log=False, tickcolor=None,
                 orientation="right", pos=[0,0], shrink=1.0):
        """
        Places a colorbar with a colormap *name* and a value range
        *zrange*, labelled with *label*.  The axis may be logged by
        setting *log*, and the tick colors are changed with
        *tickcolor*.  The *orientation* can be left, right, top, or
        bottom.  The position can be manually adjusted with *pos* and
        scaled with *shrink*.
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
            print "orientation %s unknown" % orientation
            return

        # Scale the colorbar
        shift = (0.5*(1.0-shrink)*size[0], 0.5*(1.0-shrink)*size[1])
        size = (size[0] * shrink, size[1] * shrink)
        origin = (origin[0] + shift[0], origin[1] + shift[1])

        # Convert the colormap into a string
        x = na.linspace(1,0,256)
        cm_string = cm.cmap_d[name](x, bytes=True)[:,0:3].tostring()

        cmap_im = pyx.bitmap.image(imsize[0], imsize[1], "RGB", cm_string)
        if orientation == "top" or orientation == "bottom":
            imorigin = (imorigin[0] - shift[0], imorigin[1] - shift[1])
            self.canvas.insert(pyx.bitmap.bitmap(imorigin[0], imorigin[1], cmap_im,
                                                 width=-size[0], height=size[1]))
        else:
            self.canvas.insert(pyx.bitmap.bitmap(origin[0], origin[1], cmap_im,
                                                 width=size[0], height=size[1]))

        if tickcolor is None:
            c1 = pyx.graph.axis.painter.regular\
                 (tickattrs=[pyx.color.cmyk.black])
            c2 = pyx.graph.axis.painter.regular\
                 (tickattrs=[pyx.color.cmyk.black], labelattrs=None)
        else:
            c1 = pyx.graph.axis.painter.regular(tickattrs=[tickcolor])
            c2 = pyx.graph.axis.painter.regular(tickattrs=[tickcolor],
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

    def colorbar_yt(self, plot, **kwargs):
        if plot.cmap != None:
            _cmap = plot.cmap.name
        else:
            _cmap = 'algae'
        if isinstance(plot, raven.PlotTypes.VMPlot):
            # Taken from yt
            proj = "Proj" in plot._type_name and \
                   plot.data._weight is None
            _zlabel = plot.pf.field_info[plot.axis_names["Z"]].get_label(proj)
            print _zlabel
            _zlabel = _zlabel.replace("_","\;")
            print _zlabel
            _zlog = plot.log_field
        else:
            _zlabel = plot._z_label.replace("_","\;")
            _zlog = plot._log_z
        _zrange = (plot.norm.vmin, plot.norm.vmax)
        self.colorbar(_cmap, zrange=_zrange, label=_zlabel, log=_zlog, **kwargs)

#=============================================================================

    def circle(self, radius=0.2, loc=(0.5,0.5),
               color=pyx.color.cmyk.white,
               linewidth=pyx.style.linewidth.normal):
        """
        Draws a circle with a *radius* at *loc* with a *color* and
        *linewidth*.
        """
        circle = pyx.path.circle(self.figsize[0]*loc[0],
                                 self.figsize[1]*loc[1],
                                 self.figsize[0]*radius)
        self.canvas.stroke(circle, [color, linewidth])

#=============================================================================

    def scale_line(self, size=0.2, label="", loc=(0.05,0.08), labelloc="top",
                   color=pyx.color.cmyk.white,
                   linewidth=pyx.style.linewidth.normal):
        """
        Draws a scale line with a *size*, *label* at position *loc*.
        The label position is specified with *labelloc* and can be top
        or bottom.  The *color* and *linewidth* can also be specified.
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
                  valign=pyx.text.valign.top):
        """
        Inserts a *text* box at *loc* with a text color *color* and
        background color of *bgcolor*.  The text is aligned with
        *halign* and *valign*.
        """
        tbox = self.canvas.text(self.figsize[0]*loc[0],
                                self.figsize[1]*loc[1],
                                text, [color, valign, halign])
        if bgcolor != None:
            tpath = tbox.bbox().enlarged(2*pyx.unit.x_pt).path()
            self.canvas.draw(tpath, [pyx.deco.filled([bgcolor]),
                                     pyx.deco.stroked()])
        self.canvas.insert(tbox)
        
#=============================================================================

    def save_fig(self, filename="test", format="eps"):
        """
        Saves current figure to *filename* in *format*, which can be
        eps or pdf.
        """
        if format =="eps":
            self.canvas.writeEPSfile(filename)
        elif format == "pdf":
            self.canvas.writePDFfile(filename)
        else:
            print "format %s unknown." % (format)
            
#=============================================================================
#=============================================================================
#=============================================================================

def multiplot(ncol, nrow, yt_plots=None, images=None, xranges=None,
              yranges=None, xlabels=None, ylabels=None, colorbars=None,
              shrink_cb=0.95, figsize=(8,8), margins=(0,0), titles=None,
              savefig=None, yt_nocbar=False, bare_axes=False):
    # Error check
    if images != None:
        if len(images) != ncol*nrow:
            print "Number of images (%d) doesn't match nrow(%d) x ncol(%d)." %\
                  (len(images), nrow, ncol)
            return
    if yt_plots is None and images is None:
        print "Must supply either yt_plots or image filenames."
        return
    if yt_plots != None and images != None:
        print "Given both images and yt plots.  Ignoring images."
    if yt_plots != None:
        _yt = True

    # If no ranges or labels given and given only images, fill them in.
    if not _yt:
        if xranges is None:
            xranges = []
            for i in range(nrow*ncol): xranges.append((0,1))
        if yranges is None:
            yranges = []
            for i in range(nrow*ncol): yranges.append((0,1))
        if xlabels is None:
            xlabels = []
            for i in range(nrow*ncol): xlabels.append("")
        if ylabels is None:
            ylabels = []
            for i in range(nrow*ncol): ylabels.append("")

    d = DualEPS(figsize=figsize)
    count = 0
    for j in range(nrow):
        ypos = j*(figsize[1] + margins[1])
        for i in range(ncol):
            xpos = i*(figsize[0] + margins[0])
            index = j*ncol + i
            if j == nrow-1:
                xaxis = 1
            elif j == 0:
                xaxis = 0
            else:
                xaxis = -1
            if i == 0:
                yaxis = 0
            elif i == ncol-1:
                yaxis = 1
            else:
                yaxis = -1
            if _yt:
                d.insert_image_yt(yt_plots[index], pos=(xpos, ypos))
                d.axis_box_yt(yt_plots[index], pos=(xpos, ypos),
                              bare_axes=bare_axes, xaxis_side=xaxis,
                              yaxis_side=yaxis)
            else:
                d.insert_image(images[index], pos=(xpos,ypos))
                d.axis_box(pos = (xpos, ypos),
                           xrange=xranges[index], yrange=yranges[index],
                           xlabel=xlabels[i], ylabel=ylabels[j],
                           bare_axes=bare_axes, xaxis_side=xaxis, yaxis_side=yaxis)
            if titles != None:
                if titles[index] != None:
                    d.title_box(titles[index],
                                loc=(i+0.02+i*margins[0]/figsize[0],
                                     j+0.98+j*margins[1]/figsize[1]))

    # Insert colorbars after all axes are placed because we want to
    # put them on the edges of the bounding box.
    bbox = (100.0 * d.canvas.bbox().left().t,
            100.0 * d.canvas.bbox().right().t - d.figsize[0],
            100.0 * d.canvas.bbox().bottom().t,
            100.0 * d.canvas.bbox().top().t - d.figsize[1])
    for j in range(nrow):
        ypos0 = j*(figsize[1] + margins[1])
        for i in range(ncol):
            xpos0 = i*(figsize[0] + margins[0])
            index = j*ncol + i
            if (not _yt and colorbars != None) or (_yt and not yt_nocbar):
                if _yt or colorbars[index] != None:
                    if ncol == 1:
                        orientation = "right"
                        xpos = bbox[1]
                        ypos = ypos0
                    elif i == 0:
                        orientation = "left"
                        xpos = bbox[0]
                        ypos = ypos0
                    elif i+1 == ncol:
                        orientation = "right"
                        xpos = bbox[1]
                        ypos = ypos0
                    elif j == 0:
                        orientation = "bottom"
                        ypos = bbox[2]
                        xpos = xpos0
                    elif j+1 == nrow:
                        orientation = "top"
                        ypos = bbox[3]
                        xpos = xpos0
                    else:
                        orientation = None  # Marker for interior plot

                    if orientation != None:
                        if _yt:
                            d.colorbar_yt(yt_plots[index],
                                          pos=[xpos,ypos],
                                          shrink=shrink_cb,
                                          orientation=orientation)
                        else:
                            d.colorbar(colorbars[index]["cmap"],
                                       zrange=colorbars[index]["range"],
                                       label=colorbars[index]["name"],
                                       log=colorbars[index]["log"],
                                       orientation=orientation,
                                       pos=[xpos,ypos],
                                       shrink=shrink_cb)

    if savefig != None:
        d.save_fig(savefig)

    return d

#=============================================================================

def multiplot_yt(ncol, nrow, plot_col, **kwargs):
    if len(plot_col.plots) < nrow*ncol:
        print "Number of plots in PlotCollection is less than nrow(%d) "\
              "x ncol(%d)." % (len(plot_col.plots), nrow, ncol)
        return
    figure = multiplot(ncol, nrow, yt_plots=plot_col.plots, **kwargs)
    return figure

#=============================================================================

def single_plot(plot, figsize=(12,12), cb_orient="right", bare_axes=False,
                savefig=None, file_format='eps'):
    d = DualEPS(figsize=figsize)
    d.insert_image_yt(plot)
    d.axis_box_yt(plot, bare_axes=bare_axes)
    d.colorbar_yt(plot, orientation=cb_orient)
    if savefig != None:
        d.save_fig(savefig, format=file_format)
    return d

#=============================================================================
def return_cmap(cmap="algae", label="", range=(0,1), log=False):
    return {'cmap': cmap, 'name': label, 'range': range, 'log': log}
    
#=============================================================================

#if __name__ == "__main__":
#    pf = load('/Users/jwise/runs/16Jul09_Pop3/DD0019/output_0019')
#    pc = PlotCollection(pf)
#    p = pc.add_slice('Density',0,use_colorbar=False)
#    p.set_width(0.1,'kpc')
#    p1 = pc.add_slice('Temperature',0,use_colorbar=False)
#    p1.set_width(0.1,'kpc')
#    p1.set_cmap('hot')
#    p1 = pc.add_phase_sphere(0.1, 'kpc', ['Radius', 'Density', 'H2I_Fraction'],
#                            weight='CellMassMsun')
#    p1.set_xlim(1e18,3e20)
#    p1 = pc.add_phase_sphere(0.1, 'kpc', ['Radius', 'Density', 'Temperature'],
#                            weight='CellMassMsun')
#    p1.set_xlim(1e18,3e20)
#    mp = multiplot_yt(2,2,pc,savefig="yt",shrink_cb=0.9, bare_axes=False,
#                      yt_nocbar=False, margins=(0.5,0.5))
#    #d = DualEPS()
#    #d.axis_box(xrange=(0,20.7), yrange=(0,20.7), xlabel='x [kpc]',
#    #           ylabel='y [kpc]')
#    #d.axis_box_yt(p)
#    #d.insert_image_yt(p)
#    #d.colorbar_yt(p)
#    #d.title_box("Halo 1", loc=(0.02,0.02), valign=pyx.text.valign.bottom)
#    #d.circle()
#    #d.save_fig('yt')
#    
#    images = ["density.jpg", "density.jpg", "density.jpg", "density.jpg"]
#    cbs=[]
#    cbs.append(return_cmap("algae", "Density [cm$^{-3}$]", (0,10), False))
#    cbs.append(return_cmap("jet", "HI Density", (0,5), False))
#    cbs.append(return_cmap("hot", r"Entropy [K cm$^2$]", (1e-2,1e6), True))
#    cbs.append(return_cmap("Spectral", "Stuff$_x$!", (1,300), True))
#    
#    mp = multiplot(images,2,2, margins=(0.1,0.1), titles=["1","2","3","4"],
#                   xlabels=["one","two"], ylabels=None, colorbars=cbs,
#                   shrink_cb=0.95)
#    mp.scale_line(label="$r_{vir}$", labelloc="top")
#    mp.save_fig("multiplot")

#    d = DualEPS()
#    d.axis_box(xrange=(0,250), yrange=(0,250), xlabel="x [pc]", ylabel="y [pc]",
#               xlog=False, ylog=False, tickcolor=pyx.color.cmyk.White,
#               bare_axes=False, xaxis_side=1, yaxis_side=0)
#    d.insert_image("density.jpg")
#    d.colorbar("algae", zrange=(1e-2,1e4), log=True,
#               label="Density [cm$^{-3}$]", orientation='right')
#    d.scale_line(label="$r_{vir}$", labelloc="top")
#    d.circle()
#    d.title_box(r"$z=17$")
#    d.save_fig(format='eps')

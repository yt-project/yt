"""
Callbacks to add additional functionality on to plots.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: J. S. Oishi <jsoishi@astro.berkeley.edu>
Affiliation: UC Berkeley
Author: Stephen Skory <s@skory.us>
Affiliation: UC San Diego
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk, JS Oishi, Stephen Skory.  All Rights Reserved.

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

import numpy as na

from yt.funcs import *
from _mpl_imports import *
from yt.utilities.definitions import \
    x_dict, x_names, \
    y_dict, y_names, \
    axis_names, \
    axis_labels
import _MPL

callback_registry = {}

class PlotCallback(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            callback_registry[name] = cls

    def __init__(self, *args, **kwargs):
        pass

    def convert_to_plot(self, plot, coord, offset = True):
        # coord should be a 2 x ncoord array-like datatype.
        ncoord = na.array(coord).shape[1]

        # Convert the data and plot limits to tiled numpy arrays so that
        # convert_to_plot is automatically vectorized.

        x0 = na.tile(plot.xlim[0],ncoord)
        x1 = na.tile(plot.xlim[1],ncoord)
        xx0 = na.tile(plot._axes.get_xlim()[0],ncoord)
        xx1 = na.tile(plot._axes.get_xlim()[1],ncoord)
        
        y0 = na.tile(plot.ylim[0],ncoord)
        y1 = na.tile(plot.ylim[1],ncoord)
        yy0 = na.tile(plot._axes.get_ylim()[0],ncoord)
        yy1 = na.tile(plot._axes.get_ylim()[1],ncoord)
        
        # We need a special case for when we are only given one coordinate.
        if na.array(coord).shape == (2,):
            return ((coord[0]-x0)/(x1-x0)*(xx1-xx0) + xx0,
                    (coord[1]-y0)/(y1-y0)*(yy1-yy0) + yy0)
        else:
            return ((coord[0][:]-x0)/(x1-x0)*(xx1-xx0) + xx0,
                    (coord[1][:]-y0)/(y1-y0)*(yy1-yy0) + yy0)

    def pixel_scale(self,plot):
        x0, x1 = plot.xlim
        xx0, xx1 = plot._axes.get_xlim()
        dx = (xx0 - xx1)/(x1 - x0)
        
        y0, y1 = plot.ylim
        yy0, yy1 = plot._axes.get_ylim()
        dy = (yy0 - yy1)/(y1 - y0)

        return (dx,dy)


class VelocityCallback(PlotCallback):
    _type_name = "velocity"
    def __init__(self, factor=16, scale=None, scale_units=None, normalize=False):
        """
        annotate_velocity(factor=16, scale=None, scale_units=None, normalize=False):
        
        Adds a 'quiver' plot of velocity to the plot, skipping all but
        every *factor* datapoint. *scale* is the data units per arrow
        length unit using *scale_units* (see
        matplotlib.axes.Axes.quiver for more info). if *normalize* is
        True, the velocity fields will be scaled by their local
        (in-plane) length, allowing morphological features to be more
        clearly seen for fields with substantial variation in field
        strength (normalize is not implemented and thus ignored for
        Cutting Planes).
        """
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale  = scale
        self.scale_units = scale_units
        self.normalize = normalize

    def __call__(self, plot):
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            qcb = CuttingQuiverCallback("CuttingPlaneVelocityX",
                                        "CuttingPlaneVelocityY",
                                        self.factor)
        else:
            xv = "%s-velocity" % (x_names[plot.data.axis])
            yv = "%s-velocity" % (y_names[plot.data.axis])
            qcb = QuiverCallback(xv, yv, self.factor, scale=self.scale, scale_units=self.scale_units, normalize=self.normalize)
        return qcb(plot)

class MagFieldCallback(PlotCallback):
    _type_name = "magnetic_field"
    def __init__(self, factor=16, scale=None, scale_units=None, normalize=False):
        """
        annotate_magnetic_field(factor=16, scale=None, scale_units=None, normalize=False):

        Adds a 'quiver' plot of magnetic field to the plot, skipping all but
        every *factor* datapoint. *scale* is the data units per arrow
        length unit using *scale_units* (see
        matplotlib.axes.Axes.quiver for more info). if *normalize* is
        True, the magnetic fields will be scaled by their local
        (in-plane) length, allowing morphological features to be more
        clearly seen for fields with substantial variation in field strength.
        """
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale  = scale
        self.scale_units = scale_units
        self.normalize = normalize

    def __call__(self, plot):
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            print "WARNING: Magnetic field on Cutting Plane Not implemented."
        else:
            xv = "B%s" % (x_names[plot.data.axis])
            yv = "B%s" % (y_names[plot.data.axis])
            qcb = QuiverCallback(xv, yv, self.factor, scale=self.scale, scale_units=self.scale_units, normalize=self.normalize)
        return qcb(plot)

class QuiverCallback(PlotCallback):
    _type_name = "quiver"
    def __init__(self, field_x, field_y, factor=16, scale=None, scale_units=None, normalize=False):
        """
        annotate_quiver(field_x, field_y, factor, scale=None, scale_units=None, normalize=False):

        Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
        from the associated data, skipping every *factor* datapoints
        *scale* is the data units per arrow length unit using *scale_units* 
        (see matplotlib.axes.Axes.quiver for more info)
        """
        PlotCallback.__init__(self)
        self.field_x = field_x
        self.field_y = field_y
        self.bv_x = self.bv_y = 0
        self.factor = factor
        self.scale = scale
        self.scale_units = scale_units
        self.normalize = normalize

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        nx = plot.image._A.shape[0] / self.factor
        ny = plot.image._A.shape[1] / self.factor
        pixX = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_x] - self.bv_x,
                             int(nx), int(ny),
                           (x0, x1, y0, y1),).transpose()
        pixY = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_y] - self.bv_y,
                             int(nx), int(ny),
                           (x0, x1, y0, y1),).transpose()
        X,Y = na.meshgrid(na.linspace(xx0,xx1,nx,endpoint=True),
                          na.linspace(yy0,yy1,ny,endpoint=True))
        if self.normalize:
            nn = na.sqrt(pixX**2 + pixY**2)
            pixX /= nn
            pixY /= nn
        plot._axes.quiver(X,Y, pixX, pixY, scale=self.scale, scale_units=self.scale_units)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ContourCallback(PlotCallback):
    _type_name = "contour"
    def __init__(self, field, ncont=5, factor=4, take_log=False, clim=None,
                 plot_args = None):
        """
        annotate_contour(self, field, ncont=5, factor=4, take_log=False, clim=None,
                         plot_args = None):

        Add contours in *field* to the plot.  *ncont* governs the number of
        contours generated, *factor* governs the number of points used in the
        interpolation, *take_log* governs how it is contoured and *clim* gives
        the (upper, lower) limits for contouring.
        """
        PlotCallback.__init__(self)
        self.ncont = ncont
        self.field = field
        self.factor = factor
        self.take_log = take_log
        from yt.utilities.delaunay.triangulate import Triangulation as triang
        self.triang = triang
        if self.take_log and clim is not None: clim = (na.log10(clim[0]), na.log10(clim[1]))
        if clim is not None: self.ncont = na.linspace(clim[0], clim[1], ncont)
        if plot_args is None: plot_args = {'colors':'k'}
        self.plot_args = plot_args

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        
        plot._axes.hold(True)
        
        numPoints_x = plot.image._A.shape[0]
        numPoints_y = plot.image._A.shape[1]
        
        # Multiply by dx and dy to go from data->plot
        dx = (yy1 - yy0) / (y1-y0)
        dy = (xx1 - xx0) / (x1-x0)

        #dcollins Jan 11 2009.  Improved to allow for periodic shifts in the plot.
        #Now makes a copy of the position fields "px" and "py" and adds the
        #appropriate shift to the coppied field.  

        #set the cumulative arrays for the periodic shifting.
        AllX = na.zeros(plot.data["px"].size, dtype='bool')
        AllY = na.zeros(plot.data["py"].size, dtype='bool')
        XShifted = plot.data["px"].copy()
        YShifted = plot.data["py"].copy()
        dom_x, dom_y = plot._period
        for shift in na.mgrid[-1:1:3j]:
            xlim = ((plot.data["px"] + shift*dom_x >= x0)
                 &  (plot.data["px"] + shift*dom_x <= x1))
            ylim = ((plot.data["py"] + shift*dom_y >= y0)
                 &  (plot.data["py"] + shift*dom_y <= y1))
            XShifted[xlim] += shift * dom_x
            YShifted[ylim] += shift * dom_y
            AllX |= xlim
            AllY |= ylim
        # At this point XShifted and YShifted are the shifted arrays of
        # position data in data coordinates
        wI = (AllX & AllY)

        # We want xi, yi in plot coordinates
        xi, yi = na.mgrid[xx0:xx1:numPoints_x/(self.factor*1j),\
                          yy0:yy1:numPoints_y/(self.factor*1j)]

        # This converts XShifted and YShifted into plot coordinates
        x = (XShifted[wI]-x0)*dx + xx0
        y = (YShifted[wI]-y0)*dy + yy0
        z = plot.data[self.field][wI]
        if self.take_log: z=na.log10(z)

        # Both the input and output from the triangulator are in plot
        # coordinates
        zi = self.triang(x,y).nn_interpolator(z)(xi,yi)

        plot._axes.contour(xi,yi,zi,self.ncont, **self.plot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class GridBoundaryCallback(PlotCallback):
    _type_name = "grids"
    def __init__(self, alpha=1.0, min_pix=1, annotate=False, periodic=True):
        """
        annotate_grids(alpha=1.0, min_pix=1, annotate=False, periodic=True)

        Adds grid boundaries to a plot, optionally with *alpha*-blending.
        Cuttoff for display is at *min_pix* wide.
        *annotate* puts the grid id in the corner of the grid.  (Not so great in projections...)
        """
        PlotCallback.__init__(self)
        self.alpha = alpha
        self.min_pix = min_pix
        self.annotate = annotate # put grid numbers in the corner.
        self.periodic = periodic

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        width, height = plot.image._A.shape
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        xi = x_dict[plot.data.axis]
        yi = y_dict[plot.data.axis]
        dx = width / (x1-x0)
        dy = height / (y1-y0)
        px_index = x_dict[plot.data.axis]
        py_index = y_dict[plot.data.axis]
        dom = plot.data.pf.domain_right_edge - plot.data.pf.domain_left_edge
        if self.periodic:
            pxs, pys = na.mgrid[-1:1:3j,-1:1:3j]
        else:
            pxs, pys = na.mgrid[0:0:1j,0:0:1j]
        GLE = plot.data.grid_left_edge
        GRE = plot.data.grid_right_edge
        for px_off, py_off in zip(pxs.ravel(), pys.ravel()):
            pxo = px_off * dom[px_index]
            pyo = py_off * dom[py_index]
            left_edge_px = (GLE[:,px_index]+pxo-x0)*dx
            left_edge_py = (GLE[:,py_index]+pyo-y0)*dy
            right_edge_px = (GRE[:,px_index]+pxo-x0)*dx
            right_edge_py = (GRE[:,py_index]+pyo-y0)*dy
            verts = na.array(
                [(left_edge_px, left_edge_px, right_edge_px, right_edge_px),
                 (left_edge_py, right_edge_py, right_edge_py, left_edge_py)])
            visible =  ( right_edge_px - left_edge_px > self.min_pix ) & \
                       ( right_edge_px - left_edge_px > self.min_pix )
            verts=verts.transpose()[visible,:,:]
            if verts.size == 0: continue
            edgecolors = (0.0,0.0,0.0,self.alpha)
            verts[:,:,0]= (xx1-xx0)*(verts[:,:,0]/width) + xx0
            verts[:,:,1]= (yy1-yy0)*(verts[:,:,1]/height) + yy0
            grid_collection = matplotlib.collections.PolyCollection(
                verts, facecolors="none",
                edgecolors=edgecolors)
            plot._axes.hold(True)
            plot._axes.add_collection(grid_collection)
            if self.annotate:
                ids = [g.id for g in plot.data._grids]
                for n in range(len(left_edge_px)):
                    plot._axes.text(left_edge_px[n]+2,left_edge_py[n]+2,ids[n])
            plot._axes.hold(False)

class StreamlineCallback(PlotCallback):
    _type_name = "streamlines"
    def __init__(self, field_x, field_y, factor=6.0, nx=16, ny=16,
                 xstart=(0,1), ystart=(0,1), nsample=256,
                 start_at_xedge=False, start_at_yedge=False,
                 plot_args=None):
        """
        annotate_streamlines(field_x, field_y, factor=6.0, nx=16, ny=16,
                             xstart=(0,1), ystart=(0,1), nsample=256,
                             start_at_xedge=False, start_at_yedge=False,
                             plot_args=None):

        Add streamlines to any plot, using the *field_x* and *field_y*
        from the associated data, using *nx* and *ny* starting points
        that are bounded by *xstart* and *ystart*.  To begin
        streamlines from the left edge of the plot, set
        *start_at_xedge* to True; for the bottom edge, use
        *start_at_yedge*.  A line with the qmean vector magnitude will
        cover 1.0/*factor* of the image.
        """
        PlotCallback.__init__(self)
        self.field_x = field_x
        self.field_y = field_y
        self.xstart = xstart
        self.ystart = ystart
        self.nsample = nsample
        self.factor = factor
        if start_at_xedge:
            self.data_size = (1,ny)
        elif start_at_yedge:
            self.data_size = (nx,1)
        else:
            self.data_size = (nx,ny)
        if plot_args is None: plot_args = {'color':'k', 'linestyle':'-'}
        self.plot_args = plot_args
        
    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        nx = plot.image._A.shape[0]
        ny = plot.image._A.shape[1]
        pixX = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_x],
                             int(nx), int(ny),
                           (x0, x1, y0, y1),)
        pixY = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_y],
                             int(nx), int(ny),
                           (x0, x1, y0, y1),)
        r0 = na.mgrid[self.xstart[0]*nx:self.xstart[1]*nx:self.data_size[0]*1j,
                      self.ystart[0]*ny:self.ystart[1]*ny:self.data_size[1]*1j]
        lines = na.zeros((self.nsample, 2, self.data_size[0], self.data_size[1]))
        lines[0,:,:,:] = r0
        mag = na.sqrt(pixX**2 + pixY**2)
        scale = na.sqrt(nx*ny) / (self.factor * mag.mean())
        dt = 1.0 / (self.nsample-1)
        for i in range(1,self.nsample):
            xt = lines[i-1,0,:,:]
            yt = lines[i-1,1,:,:]
            ix = na.maximum(na.minimum((xt).astype('int'), nx-1), 0)
            iy = na.maximum(na.minimum((yt).astype('int'), ny-1), 0)
            lines[i,0,:,:] = xt + dt * pixX[ix,iy] * scale
            lines[i,1,:,:] = yt + dt * pixY[ix,iy] * scale
        for i in range(self.data_size[0]):
            for j in range(self.data_size[1]):
                plot._axes.plot(lines[:,0,i,j], lines[:,1,i,j],
                                **self.plot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class LabelCallback(PlotCallback):
    _type_name = "axis_label"
    def __init__(self, label):
        """
        This adds a label to the plot.
        """
        PlotCallback.__init__(self)
        self.label = label

    def __call__(self, plot):
        plot._figure.subplots_adjust(hspace=0, wspace=0,
                                     bottom=0.1, top=0.9,
                                     left=0.0, right=1.0)
        plot._axes.set_xlabel(self.label)
        plot._axes.set_ylabel(self.label)

def get_smallest_appropriate_unit(v, pf):
    max_nu = 1e30
    good_u = None
    for unit in ['mpc','kpc','pc','au','rsun','cm']:
        vv = v*pf[unit]
        if vv < max_nu and vv > 1.0:
            good_u = unit
            max_nu = v*pf[unit]
    if good_u is None : good_u = 'cm'
    return good_u

class UnitBoundaryCallback(PlotCallback):
    _type_name = "units"
    def __init__(self, unit = "au", factor=4, text_annotate=True, text_which=-2):
        """
        Add on a plot indicating where *factor*s of *unit* are shown.
        Optionally *text_annotate* on the *text_which*-indexed box on display.
        """
        PlotCallback.__init__(self)
        self.unit = unit
        self.factor = factor
        self.text_annotate = text_annotate
        self.text_which = -2

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = mpl_get_bounds(plot._axes.bbox)
        xi = x_dict[plot.data.axis]
        yi = y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        center = plot.data.center
        min_dx = plot.data['pdx'].min()
        max_dx = plot.data['pdx'].max()
        w_min_x = 250.0 * min_dx
        w_max_x = 1.0 / self.factor
        min_exp_x = na.ceil(na.log10(w_min_x*plot.data.pf[self.unit])
                           /na.log10(self.factor))
        max_exp_x = na.floor(na.log10(w_max_x*plot.data.pf[self.unit])
                            /na.log10(self.factor))
        n_x = max_exp_x - min_exp_x + 1
        widths = na.logspace(min_exp_x, max_exp_x, num = n_x, base=self.factor)
        widths /= plot.data.pf[self.unit]
        left_edge_px = (center[xi] - widths/2.0 - x0)*dx
        left_edge_py = (center[yi] - widths/2.0 - y0)*dy
        right_edge_px = (center[xi] + widths/2.0 - x0)*dx
        right_edge_py = (center[yi] + widths/2.0 - y0)*dy
        verts = na.array(
                [(left_edge_px, left_edge_px, right_edge_px, right_edge_px),
                 (left_edge_py, right_edge_py, right_edge_py, left_edge_py)])
        visible =  ( right_edge_px - left_edge_px > 25 ) & \
                   ( right_edge_px - left_edge_px > 25 ) & \
                   ( (right_edge_px < width) & (left_edge_px > 0) ) & \
                   ( (right_edge_py < height) & (left_edge_py > 0) )
        verts=verts.transpose()[visible,:,:]
        grid_collection = matplotlib.collections.PolyCollection(
                verts, facecolors="none",
                       edgecolors=(0.0,0.0,0.0,1.0),
                       linewidths=2.5)
        plot._axes.hold(True)
        plot._axes.add_collection(grid_collection)
        if self.text_annotate:
            ti = max(self.text_which, -1*len(widths[visible]))
            if ti < len(widths[visible]): 
                w = widths[visible][ti]
                good_u = get_smallest_appropriate_unit(w, plot.data.pf)
                w *= plot.data.pf[good_u]
                plot._axes.annotate("%0.3e %s" % (w,good_u), verts[ti,1,:]+5)
        plot._axes.hold(False)

class LinePlotCallback(PlotCallback):
    _type_name = "line"
    def __init__(self, x, y, plot_args = None):
        """
        annotate_line(x, y, plot_args = None)

        Over plot *x* and *y* with *plot_args* fed into the plot.
        """
        PlotCallback.__init__(self)
        self.x = x
        self.y = y
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        plot._axes.hold(True)
        plot._axes.plot(self.x, self.y, **self.plot_args)
        plot._axes.hold(False)

class ImageLineCallback(LinePlotCallback):
    _type_name = "image_line"

    def __init__(self, p1, p2, data_coords=False, plot_args = None):
        """
        annotate_image_line(p1, p2, data_coords=False, plot_args = None)

        Plot from *p1* to *p2* (image plane coordinates)
        with *plot_args* fed into the plot.
        """
        PlotCallback.__init__(self)
        self.p1 = p1
        self.p2 = p2
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args
        self._ids = []
        self.data_coords = data_coords

    def __call__(self, plot):
        # We manually clear out any previous calls to this callback:
        plot._axes.lines = [l for l in plot._axes.lines if id(l) not in self._ids]
        kwargs = self.plot_args.copy()
        if self.data_coords and len(plot.image._A.shape) == 2:
            p1 = self.convert_to_plot(plot, self.p1)
            p2 = self.convert_to_plot(plot, self.p2)
        else:
            p1, p2 = self.p1, self.p2
            if not self.data_coords:
                kwargs["transform"] = plot._axes.transAxes

        px, py = (p1[0], p2[0]), (p1[1], p2[1])

        # Save state
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        ii = plot._axes.plot(px, py, **kwargs)
        self._ids.append(id(ii[0]))
        # Reset state
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class CuttingQuiverCallback(PlotCallback):
    _type_name = "cquiver"
    def __init__(self, field_x, field_y, factor):
        """
        annotate_cquiver(field_x, field_y, factor)

        Get a quiver plot on top of a cutting plane, using *field_x* and
        *field_y*, skipping every *factor* datapoint in the discretization.
        """
        PlotCallback.__init__(self)
        self.field_x = field_x
        self.field_y = field_y
        self.factor = factor

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        nx = plot.image._A.shape[0] / self.factor
        ny = plot.image._A.shape[1] / self.factor
        indices = na.argsort(plot.data['dx'])[::-1]
        pixX = _MPL.CPixelize( plot.data['x'], plot.data['y'], plot.data['z'],
                               plot.data['px'], plot.data['py'],
                               plot.data['pdx'], plot.data['pdy'], plot.data['pdz'],
                               plot.data.center, plot.data._inv_mat, indices,
                               plot.data[self.field_x],
                               int(nx), int(ny),
                               (x0, x1, y0, y1),).transpose()
        pixY = _MPL.CPixelize( plot.data['x'], plot.data['y'], plot.data['z'],
                               plot.data['px'], plot.data['py'],
                               plot.data['pdx'], plot.data['pdy'], plot.data['pdz'],
                               plot.data.center, plot.data._inv_mat, indices,
                               plot.data[self.field_y],
                               int(nx), int(ny),
                               (x0, x1, y0, y1),).transpose()
        X = na.mgrid[0:plot.image._A.shape[0]-1:nx*1j]# + 0.5*factor
        Y = na.mgrid[0:plot.image._A.shape[1]-1:ny*1j]# + 0.5*factor
        plot._axes.quiver(X,Y, pixX, pixY)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ClumpContourCallback(PlotCallback):
    _type_name = "clumps"
    def __init__(self, clumps, plot_args = None):
        """
        annotate_clumps(clumps, plot_args = None)

        Take a list of *clumps* and plot them as a set of contours.
        """
        self.clumps = clumps
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)

        px_index = x_dict[plot.data.axis]
        py_index = y_dict[plot.data.axis]

        xf = axis_names[px_index]
        yf = axis_names[py_index]

        DomainRight = plot.data.pf.domain_right_edge
        DomainLeft = plot.data.pf.domain_left_edge
        DomainWidth = DomainRight - DomainLeft
        
        nx, ny = plot.image._A.shape
        buff = na.zeros((nx,ny),dtype='float64')
        for i,clump in enumerate(reversed(self.clumps)):
            mylog.debug("Pixelizing contour %s", i)


            xf_copy = clump[xf].copy()
            yf_copy = clump[yf].copy()
            
            temp = _MPL.Pixelize(xf_copy, yf_copy, 
                                 clump['dx']/2.0,
                                 clump['dy']/2.0,
                                 clump['dx']*0.0+i+1, # inits inside Pixelize
                                 int(nx), int(ny),
                             (x0, x1, y0, y1), 0).transpose()
            buff = na.maximum(temp, buff)
        self.rv = plot._axes.contour(buff, len(self.clumps)+1,
                                     **self.plot_args)
        plot._axes.hold(False)

class ArrowCallback(PlotCallback):
    _type_name = "arrow"
    def __init__(self, pos, code_size, plot_args = None):
        """
        annotate_arrow(pos, code_size, plot_args = None)

        This adds an arrow pointing at *pos* with size *code_size* in code
        units.  *plot_args* is a dict fed to matplotlib with arrow properties.
        """
        self.pos = pos
        if not iterable(code_size):
            code_size = (code_size, code_size)
        self.code_size = code_size
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        from matplotlib.patches import Arrow
        # Now convert the pixels to code information
        x, y = self.convert_to_plot(plot, self.pos)
        dx, dy = self.convert_to_plot(plot, self.code_size, False)
        arrow = Arrow(x, y, dx, dy, **self.plot_args)
        plot._axes.add_patch(arrow)

class PointAnnotateCallback(PlotCallback):
    _type_name = "point"
    def __init__(self, pos, text, text_args = None):
        """
        annotate_point(pos, text, text_args = None)

        This adds *text* at position *pos*, where *pos* is in code-space.
        *text_args* is a dict fed to the text placement code.
        """
        self.pos = pos
        self.text = text
        if text_args is None: text_args = {}
        self.text_args = text_args

    def __call__(self, plot):


        width,height = plot.image._A.shape
        x,y = self.convert_to_plot(plot, self.pos)
        x,y = x/width,y/height

        plot._axes.text(x, y, self.text, **self.text_args)

class MarkerAnnotateCallback(PlotCallback):
    _type_name = "marker"
    def __init__(self, pos, marker='x', plot_args=None):
        """
        annotate_marker(pos, marker='x', plot_args=None)

        Adds text *marker* at *pos* in code-arguments.  *plot_args* is a dict
        that will be forwarded to the plot command.
        """
        self.pos = pos
        self.marker = marker
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        if len(self.pos) == 3:
            pos = (self.pos[x_dict[plot.data.axis]],
                   self.pos[y_dict[plot.data.axis]])
        else: pos = self.pos
        x,y = self.convert_to_plot(plot, pos)
        plot._axes.hold(True)
        plot._axes.plot((x,),(y,),self.marker, **self.plot_args)
        plot._axes.hold(False)

class SphereCallback(PlotCallback):
    _type_name = "sphere"
    def __init__(self, center, radius, circle_args = None,
                 text = None, text_args = None):
        """
        annotate_sphere(center, radius, circle_args = None,
                        text = None, text_args = None)
        
        A sphere centered at *center* in code units with radius *radius* in
        code units will be created, with optional *circle_args*, *text*, and
        *text_args*.
        """
        self.center = center
        self.radius = radius
        if circle_args is None: circle_args = {}
        if 'fill' not in circle_args: circle_args['fill'] = False
        self.circle_args = circle_args
        self.text = text
        self.text_args = text_args
        if self.text_args is None: self.text_args = {}

    def __call__(self, plot):
        from matplotlib.patches import Circle
        
        radius = self.radius * self.pixel_scale(plot)[0]

        (xi, yi) = (x_dict[plot.data.axis], y_dict[plot.data.axis])

        (center_x,center_y) = self.convert_to_plot(plot,(self.center[xi], self.center[yi]))
        
        cir = Circle((center_x, center_y), radius, **self.circle_args)
        plot._axes.add_patch(cir)
        if self.text is not None:
            plot._axes.text(center_x, center_y, self.text,
                            **self.text_args)

class HopCircleCallback(PlotCallback):
    _type_name = "hop_circles"
    def __init__(self, hop_output, max_number=None,
                 annotate=False, min_size=20, max_size=10000000,
                 font_size=8, print_halo_size=False,
                 print_halo_mass=False, width=None):
        """
        annotate_hop_circles(hop_output, max_number=None,
                             annotate=False, min_size=20, max_size=10000000,
                             font_size=8, print_halo_size=False,
                             print_halo_mass=False, width=None)

        Accepts a :class:`yt.HopList` *hop_output* and plots up to
        *max_number* (None for unlimited) halos as circles.
        """
        self.hop_output = hop_output
        self.max_number = max_number
        self.annotate = annotate
        self.min_size = min_size
        self.max_size = max_size
        self.font_size = font_size
        self.print_halo_size = print_halo_size
        self.print_halo_mass = print_halo_mass
        self.width = width

    def __call__(self, plot):
        from matplotlib.patches import Circle
        for halo in self.hop_output[:self.max_number]:
            size = halo.get_size()
            if size < self.min_size or size > self.max_size: continue
            # This could use halo.maximum_radius() instead of width
            if self.width is not None and \
                na.abs(halo.center_of_mass() - 
                       plot.data.center)[plot.data.axis] > \
                   self.width:
                continue
            
            radius = halo.maximum_radius() * self.pixel_scale(plot)[0]
            center = halo.center_of_mass()
            
            (xi, yi) = (x_dict[plot.data.axis], y_dict[plot.data.axis])

            (center_x,center_y) = self.convert_to_plot(plot,(center[xi], center[yi]))
            cir = Circle((center_x, center_y), radius, fill=False)
            plot._axes.add_patch(cir)
            if self.annotate:
                if self.print_halo_size:
                    plot._axes.text(center_x, center_y, "%s" % size,
                    fontsize=self.font_size)
                elif self.print_halo_mass:
                    plot._axes.text(center_x, center_y, "%s" % halo.total_mass(),
                    fontsize=self.font_size)
                else:
                    plot._axes.text(center_x, center_y, "%s" % halo.id,
                    fontsize=self.font_size)

class HopParticleCallback(PlotCallback):
    _type_name = "hop_particles"
    def __init__(self, hop_output, max_number=None, p_size=1.0,
                min_size=20, alpha=0.2):
        """
        annotate_hop_particles(hop_output, max_number, p_size=1.0,
                               min_size=20, alpha=0.2):

        Adds particle positions for the members of each halo as identified
        by HOP. Along *axis* up to *max_number* groups in *hop_output* that are
        larger than *min_size* are plotted with *p_size* pixels per particle; 
        *alpha* determines the opacity of each particle.
        """
        self.hop_output = hop_output
        self.p_size = p_size
        self.max_number = max_number
        self.min_size = min_size
        self.alpha = alpha
    
    def __call__(self,plot):
        (dx,dy) = self.pixel_scale(plot)

        (xi, yi) = (x_names[plot.data.axis], y_names[plot.data.axis])

        # now we loop over the haloes
        for halo in self.hop_output[:self.max_number]:
            size = halo.get_size()

            if size < self.min_size: continue

            (px,py) = self.convert_to_plot(plot,(halo["particle_position_%s" % xi],
                                                 halo["particle_position_%s" % yi]))
            
            # Need to get the plot limits and set the hold state before scatter
            # and then restore the limits and turn off the hold state afterwards
            # because scatter will automatically adjust the plot window which we
            # do not want
            
            xlim = plot._axes.get_xlim()
            ylim = plot._axes.get_ylim()
            plot._axes.hold(True)

            plot._axes.scatter(px, py, edgecolors="None",
                s=self.p_size, c='black', alpha=self.alpha)
            
            plot._axes.set_xlim(xlim)
            plot._axes.set_ylim(ylim)
            plot._axes.hold(False)


class CoordAxesCallback(PlotCallback):
    _type_name = "coord_axes"
    def __init__(self,unit=None,coords=False):
        """
        Creates x and y axes for a VMPlot. In the future, it will
        attempt to guess the proper units to use.
        """
        PlotCallback.__init__(self)
        self.unit = unit
        self.coords = coords

    def __call__(self,plot):
        # 1. find out what the domain is
        # 2. pick a unit for it
        # 3. run self._axes.set_xlabel & self._axes.set_ylabel to actually lay things down.
        # 4. adjust extent information to make sure labels are visable.

        # put the plot into data coordinates
        nx,ny = plot.image._A.shape
        dx = (plot.xlim[1] - plot.xlim[0])/nx
        dy = (plot.ylim[1] - plot.ylim[0])/ny

        unit_conversion = plot.pf[plot.im["Unit"]]
        aspect = (plot.xlim[1]-plot.xlim[0])/(plot.ylim[1]-plot.ylim[0])

        print "aspect ratio = ", aspect

        # if coords is False, label axes relative to the center of the
        # display. if coords is True, label axes with the absolute
        # coordinates of the region.
        xcenter = 0.
        ycenter = 0.
        if not self.coords:
            center = plot.data.center
            if plot.data.axis == 0:
                xcenter = center[1]
                ycenter = center[2]
            elif plot.data.axis == 1:
                xcenter = center[0]
                ycenter = center[2]
            else:
                xcenter = center[0]
                ycenter = center[1]


            xformat_function = lambda a,b: '%7.1e' %((a*dx + plot.xlim[0] - xcenter)*unit_conversion)
            yformat_function = lambda a,b: '%7.1e' %((a*dy + plot.ylim[0] - ycenter)*unit_conversion)
        else:
            xformat_function = lambda a,b: '%7.1e' %((a*dx + plot.xlim[0])*unit_conversion)
            yformat_function = lambda a,b: '%7.1e' %((a*dy + plot.ylim[0])*unit_conversion)
            
        xticker = matplotlib.ticker.FuncFormatter(xformat_function)
        yticker = matplotlib.ticker.FuncFormatter(yformat_function)
        plot._axes.xaxis.set_major_formatter(xticker)
        plot._axes.yaxis.set_major_formatter(yticker)
        
        xlabel = '%s (%s)' % (axis_labels[plot.data.axis][0],plot.im["Unit"])
        ylabel = '%s (%s)' % (axis_labels[plot.data.axis][1],plot.im["Unit"])
        xticksize = nx/4.
        yticksize = ny/4.
        plot._axes.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([i*xticksize for i in range(0,5)]))
        plot._axes.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([i*yticksize for i in range(0,5)]))
        
        plot._axes.set_xlabel(xlabel,visible=True)
        plot._axes.set_ylabel(ylabel,visible=True)
        plot._figure.subplots_adjust(left=0.1,right=0.8)

class TextLabelCallback(PlotCallback):
    _type_name = "text"
    def __init__(self, pos, text, data_coords=False, text_args = None):
        """
        annotate_text(pos, text, data_coords=False, text_args = None)

        Accepts a position in (0..1, 0..1) of the image, some text and
        optionally some text arguments. If data_coords is True,
        position will be in code units instead of image coordinates.
        """
        self.pos = pos
        self.text = text
        self.data_coords = data_coords
        if text_args is None: text_args = {}
        self.text_args = text_args

    def __call__(self, plot):
        kwargs = self.text_args.copy()
        if self.data_coords and len(plot.image._A.shape) == 2:
            if len(self.pos) == 3:
                pos = (self.pos[x_dict[plot.data.axis]],
                       self.pos[y_dict[plot.data.axis]])
            else: pos = self.pos
            x,y = self.convert_to_plot(plot, pos)
        else:
            x, y = self.pos
            if not self.data_coords:
                kwargs["transform"] = plot._axes.transAxes
        plot._axes.text(x, y, self.text, **kwargs)

class ParticleCallback(PlotCallback):
    _type_name = "particles"
    region = None
    _descriptor = None
    def __init__(self, width, p_size=1.0, col='k', marker='o', stride=1.0,
                 ptype=None, stars_only=False, dm_only=False,
                 minimum_mass=None, alpha=1.0):
        """
        annotate_particles(width, p_size=1.0, col='k', marker='o', stride=1.0,
                           ptype=None, stars_only=False, dm_only=False,
                           minimum_mass=None, alpha=1.0)

        Adds particle positions, based on a thick slab along *axis* with a
        *width* along the line of sight.  *p_size* controls the number of
        pixels per particle, and *col* governs the color.  *ptype* will
        restrict plotted particles to only those that are of a given type.
        *minimum_mass* will require that the particles be of a given mass,
        calculated via ParticleMassMsun, to be plotted. *alpha* determines
        each particle's opacity.
        """
        PlotCallback.__init__(self)
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.stars_only = stars_only
        self.dm_only = dm_only
        self.minimum_mass = minimum_mass
        self.alpha = alpha

    def __call__(self, plot):
        data = plot.data
        # we construct a recantangular prism
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        reg = self._get_region((x0,x1), (y0,y1), plot.data.axis, data)
        field_x = "particle_position_%s" % axis_names[x_dict[data.axis]]
        field_y = "particle_position_%s" % axis_names[y_dict[data.axis]]
        gg = ( ( reg[field_x] >= x0 ) & ( reg[field_x] <= x1 )
           &   ( reg[field_y] >= y0 ) & ( reg[field_y] <= y1 ) )
        if self.ptype is not None:
            gg &= (reg["particle_type"] == self.ptype)
            if gg.sum() == 0: return
        if self.stars_only:
            gg &= (reg["creation_time"] > 0.0)
            if gg.sum() == 0: return
        if self.dm_only:
            gg &= (reg["creation_time"] <= 0.0)
            if gg.sum() == 0: return
        if self.minimum_mass is not None:
            gg &= (reg["ParticleMassMsun"] >= self.minimum_mass)
            if gg.sum() == 0: return
        plot._axes.hold(True)
        px, py = self.convert_to_plot(plot,
                    [reg[field_x][gg][::self.stride],
                     reg[field_y][gg][::self.stride]])
        plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
                           s=self.p_size, c=self.color,alpha=self.alpha)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)


    def _get_region(self, xlim, ylim, axis, data):
        LE, RE = [None]*3, [None]*3
        xax = x_dict[axis]
        yax = y_dict[axis]
        zax = axis
        LE[xax], RE[xax] = xlim
        LE[yax], RE[yax] = ylim
        LE[zax] = data.center[zax] - self.width*0.5
        RE[zax] = data.center[zax] + self.width*0.5
        if self.region is not None \
            and na.all(self.region.left_edge <= LE) \
            and na.all(self.region.right_edge >= RE):
            return self.region
        self.region = data.pf.h.periodic_region(
            data.center, LE, RE)
        return self.region

class TitleCallback(PlotCallback):
    _type_name = "title"
    def __init__(self, title):
        """
        annotate_title(title)

        Accepts a *title* and adds it to the plot
        """
        PlotCallback.__init__(self)
        self.title = title

    def __call__(self,plot):
        plot._axes.set_title(self.title)


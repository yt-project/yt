"""
Callbacks to add additional functionality on to plots.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: J. S. Oishi <jsoishi@astro.berkeley.edu>
Affiliation: UC Berkeley
Author: Stephen Skory <sskory@physics.ucsd.edu>
Affiliation: UC San Diego
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk, JS Oishi, Stephen Skory.  All Rights Reserved.

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

from yt.raven import *
from PlotTypes import _get_bounds

import _MPL
import copy
callback_registry = []

class PlotCallback(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            callback_registry.append((name, cls))

    def __init__(self, *args, **kwargs):
        pass

    def convert_to_pixels(self, plot, coord, offset = True):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = _get_bounds(plot._axes.bbox)
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        return ((coord[0] - int(offset)*x0)*dx,
                (coord[1] - int(offset)*y0)*dy)

class QuiverCallback(PlotCallback):
    def __init__(self, field_x, field_y, factor):
        """
        Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
        from the associated data, skipping every *factor* datapoints.
        """
        PlotCallback.__init__(self)
        self.field_x = field_x
        self.field_y = field_y
        self.bv_x = self.bv_y = 0
        self.factor = factor

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
        X = na.mgrid[0:plot.image._A.shape[0]-1:nx*1j]# + 0.5*factor
        Y = na.mgrid[0:plot.image._A.shape[1]-1:ny*1j]# + 0.5*factor
        plot._axes.quiver(X,Y, pixX, pixY)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ParticleCallback(PlotCallback):
    def __init__(self, axis, width, p_size=1.0, col='k', stride=1.0):
        """
        Adds particle positions, based on a thick slab along *axis* with a
        *width* along the line of sight.  *p_size* controls the number of
        pixels per particle, and *col* governs the color.
        """
        PlotCallback.__init__(self)
        self.axis = axis
        self.width = width
        self.p_size = p_size
        self.color = col
        if check_color(col):
            self.color_field = False
        else:
            self.color_field = True
        self.field_x = "particle_position_%s" % lagos.axis_names[lagos.x_dict[axis]]
        self.field_y = "particle_position_%s" % lagos.axis_names[lagos.y_dict[axis]]
        self.field_z = "particle_position_%s" % lagos.axis_names[axis]
        self.stride = stride
        self.particles_x = self.particles_y = self.particles_z = None

    def _setup_particles(self, data):
        if self.particles_x is not None: return
        grids = data._grids
        particles_x = []; particles_y = []; particles_z = [];
        for g in grids:
            particles_x.append(g[self.field_x][::self.stride])
            particles_y.append(g[self.field_y][::self.stride])
            particles_z.append(g[self.field_z][::self.stride])
        self.particles_x = na.concatenate(particles_x); del particles_x
        self.particles_y = na.concatenate(particles_y); del particles_y
        self.particles_z = na.concatenate(particles_z); del particles_z
        if not self.color_field: return
        particles_c = []
        for g in grids:
            particles_c.append(g[self.color][::self.stride])
        self.particles_c = na.log10(na.concatenate(particles_c)); del particles_c

    def __call__(self, plot):
        z0 = plot.data.center[self.axis] - self.width/2.0
        z1 = plot.data.center[self.axis] + self.width/2.0
        self._setup_particles(plot.data)
        if len(self.particles_x) == 0: return
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        print "Particle bounding box:", x0, x1, y0, y1, z0, z1
        # Now we rescale because our axes limits != data limits
        goodI = na.where( (self.particles_x < x1) & (self.particles_x > x0)
                        & (self.particles_y < y1) & (self.particles_y > y0)
                        & (self.particles_z < z1) & (self.particles_z > z0))
        particles_x = (self.particles_x[goodI] - x0) * (xx1-xx0)/(x1-x0) + xx0
        particles_y = (self.particles_y[goodI] - y0) * (yy1-yy0)/(y1-y0) + yy0
        print "Particle px extrema", particles_x.min(), particles_x.max(), \
                                     particles_y.min(), particles_y.max()
        print "Axial limits", xx0, xx1, yy0, yy1
        if not self.color_field: particles_c = self.color
        else: particles_c = self.particles_c[goodI]
        plot._axes.hold(True)
        plot._axes.scatter(particles_x, particles_y, edgecolors='None',
                           s=self.p_size, c=particles_c)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ContourCallback(PlotCallback):
    def __init__(self, field, ncont=5, factor=4, take_log=False, clim=None,
                 plot_args = None):
        """
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
        try:
            import delaunay as de
            self.de = de
        except ImportError:
            mylog.warning("Callback failed; no delaunay module")
            self.__call__ = lambda a: None
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
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        #dcollins Jan 11 2009.  Improved to allow for periodic shifts in the plot.
        #Now makes a copy of the position fields "px" and "py" and adds the
        #appropriate shift to the coppied field.  
        DomainWidth = plot.data.pf["DomainRightEdge"] - plot.data.pf["DomainLeftEdge"]
        px_index = lagos.x_dict[plot.data.axis]
        py_index = lagos.y_dict[plot.data.axis]

        #set the cumulative arrays for the periodic shifting.
        AllX = na.array([False]*plot.data["px"].size)
        AllY = na.array([False]*plot.data["py"].size)
        XShifted = copy.copy(plot.data["px"])
        YShifted = copy.copy(plot.data["py"])
        for shift in na.mgrid[-1:1:3j]*DomainWidth:
            xlim = na.logical_and(plot.data["px"] + shift >= x0*0.9,
                                  plot.data["px"] + shift <= x1*1.1)
            ylim = na.logical_and(plot.data["py"] + shift >= y0*0.9,
                                  plot.data["py"] + shift <= y1*1.1)

            XShifted[na.where(xlim)] += shift
            YShifted[na.where(ylim)] += shift
            AllX = na.logical_or(AllX, xlim)
            AllY = na.logical_or(AllY, ylim)
        wI = na.where(na.logical_and(AllX,AllY))
        xi, yi = na.mgrid[0:numPoints_x:numPoints_x/(self.factor*1j),\
                          0:numPoints_y:numPoints_y/(self.factor*1j)]
        x = (XShifted[wI]-x0)*dx 
        y = (YShifted[wI]-y0)*dy
        z = plot.data[self.field][wI]
        if self.take_log: z=na.log10(z)
        zi = self.de.Triangulation(x,y).nn_interpolator(z)(xi,yi)
        print z.min(), z.max(), na.nanmin(z), na.nanmax(z)
        print zi.min(), zi.max(), na.nanmin(zi), na.nanmax(zi)
        plot._axes.contour(xi,yi,zi,self.ncont, **self.plot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class GridBoundaryCallback(PlotCallback):
    def __init__(self, alpha=1.0, min_pix = 1):
        """
        Adds grid boundaries to a plot, optionally with *alpha*-blending.
        Cuttoff for display is at *min_pix* wide.
        """
        PlotCallback.__init__(self)
        self.alpha = alpha
        self.min_pix = min_pix

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        dx = (xx1-xx0)/(x1-x0)
        dy = (yy1-yy0)/(y1-y0)
        px_index = lagos.x_dict[plot.data.axis]
        py_index = lagos.y_dict[plot.data.axis]
        dom = plot.data.pf["DomainRightEdge"] - plot.data.pf["DomainLeftEdge"]
        for px_off, py_off in na.mgrid[-1:1:3j,-1:1:3j]:
            GLE = plot.data.gridLeftEdge + px_off * dom[px_index]
            GRE = plot.data.gridRightEdge + py_off * dom[py_index]
            left_edge_px = na.maximum((GLE[:,px_index]-x0)*dx, xx0)
            left_edge_py = na.maximum((GLE[:,py_index]-y0)*dy, yy0)
            right_edge_px = na.minimum((GRE[:,px_index]-x0)*dx, xx1)
            right_edge_py = na.minimum((GRE[:,py_index]-y0)*dy, yy1)
            verts = na.array(
                    [(left_edge_px, left_edge_px, right_edge_px, right_edge_px),
                     (left_edge_py, right_edge_py, right_edge_py, left_edge_py)])
            visible =  ( right_edge_px - left_edge_px > self.min_pix ) & \
                       ( right_edge_px - left_edge_px > self.min_pix )
            verts=verts.transpose()[visible,:,:]
            if verts.size == 0: continue
            edgecolors = (0.0,0.0,0.0,self.alpha)
            grid_collection = matplotlib.collections.PolyCollection(
                    verts, facecolors=(0.0,0.0,0.0,0.0),
                           edgecolors=edgecolors)
            plot._axes.hold(True)
            plot._axes.add_collection(grid_collection)
            plot._axes.hold(False)

class LabelCallback(PlotCallback):
    def __init__(self, label):
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
    return good_u

class UnitBoundaryCallback(PlotCallback):
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
        l, b, width, height = _get_bounds(plot._axes.bbox)
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
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
                verts, facecolors=(0.0,0.0,0.0,0.0),
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
    def __init__(self, x, y, plot_args = None):
        """
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

class CuttingQuiverCallback(PlotCallback):
    def __init__(self, field_x, field_y, factor):
        """
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
    def __init__(self, clumps, axis = None, plot_args = None):
        """
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

        px_index = lagos.x_dict[plot.data.axis]
        py_index = lagos.y_dict[plot.data.axis]

        xf = lagos.axis_names[px_index]
        yf = lagos.axis_names[py_index]

        DomainRight = plot.data.pf["DomainRightEdge"]
        DomainLeft = plot.data.pf["DomainLeftEdge"]
        DomainWidth = DomainRight - DomainLeft
        
        nx, ny = plot.image._A.shape
        buff = na.zeros((nx,ny),dtype='float64')
        for i,clump in enumerate(reversed(self.clumps)):
            mylog.debug("Pixelizing contour %s", i)


            xf_copy = copy.copy(clump[xf])
            yf_copy = copy.copy(clump[yf])

            #Shift zones that belong shifted, both directions in X and Y.
            shifted = na.logical_and( xf_copy + DomainWidth[px_index] >= DomainRight[px_index],
                                      xf_copy + DomainWidth[px_index]<= x1 )
            xf_copy[na.where(shifted)] += DomainWidth[px_index]
            
            shifted = na.logical_and( xf_copy - DomainWidth[px_index] <= DomainLeft[px_index],
                                      xf_copy - DomainWidth[px_index] >= x0 )
            xf_copy[na.where(shifted)] -= DomainWidth[px_index]
            
            shifted = na.logical_and( yf_copy + DomainWidth[py_index] >= DomainRight[py_index],
                                      yf_copy + DomainWidth[py_index] <= y1 )
            yf_copy[na.where(shifted)] += DomainWidth[py_index]
            
            shifted = na.logical_and( yf_copy - DomainWidth[py_index] <= DomainLeft[py_index],
                                      yf_copy - DomainWidth[py_index] >= y0 )
            yf_copy[na.where(shifted)] -= DomainWidth[py_index]
            
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
    def __init__(self, pos, code_size, plot_args = None):
        self.pos = pos
        self.code_size = code_size
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        from matplotlib.patches import Arrow
        # Now convert the pixels to code information
        x, y = self.convert_to_pixels(plot, self.pos)
        dx, dy = self.convert_to_pixels(plot, self.code_size, False)
        arrow = Arrow(x, y, dx, dy, **self.plot_args)
        plot._axes.add_patch(arrow)

class PointAnnotateCallback(PlotCallback):
    def __init__(self, pos, text, text_args = None):
        self.pos = pos
        self.text = text
        self.text_args = text_args

    def __call__(self, plot):
        x,y = self.convert_to_pixels(plot, self.pos)
        plot._axes.text(x, y, self.text, **self.text_args)

class SphereCallback(PlotCallback):
    def __init__(self, center, radius, circle_args = None,
                 text = None, text_args = None):
        self.center = center
        self.radius = radius
        self.circle_args = circle_args
        self.text = text
        self.text_args = text_args

    def __call__(self, plot):
        from matplotlib.patches import Circle
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = _get_bounds(plot._axes.bbox)
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        radius = self.radius * dx
        center_x = (self.center[xi] - x0)*dx
        center_y = (self.center[yi] - y0)*dy
        # origin = lower?  not sure why center_y and center_x are reversed
        cir = Circle((center_y, center_x), radius, fill=False,
                     **self.circle_args)
        plot._axes.add_patch(cir)
        if self.text is not None:
            plot._axes.text(center_x, center_y, "%s" % halo.id,
                            **self.text_args)

        

class HopCircleCallback(PlotCallback):
    def __init__(self, hop_output, axis, max_number=None,
                 annotate=False, min_size=20, font_size=8, print_halo_size=False,
                 print_halo_mass=False):
        self.axis = axis
        self.hop_output = hop_output
        self.max_number = max_number
        self.annotate = annotate
        self.min_size = min_size
        self.font_size = font_size
        self.print_halo_size = print_halo_size
        self.print_halo_mass = print_halo_mass

    def __call__(self, plot):
        from matplotlib.patches import Circle
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = _get_bounds(plot._axes.bbox)
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        for halo in self.hop_output[:self.max_number]:
            size = halo.get_size()
            if size < self.min_size: continue
            radius = halo.maximum_radius() * dx
            center = halo.center_of_mass()
            center_x = (center[xi] - x0)*dx
            center_y = (center[yi] - y0)*dy
            cir = Circle((center_x, center_y), radius, fill=False)
            plot._axes.add_patch(cir)
            if self.annotate:
                if self.print_halo_size:
                    plot._axes.text(center_x, center_y, "%s" % size,
                    fontsize=self.font_size)
                elif self.print_halo_mass:
                    plot._axes_text(center_x, center_y, "%s" % halo.total_mass(),
                    fontsize=self.font_size)
                else:
                    plot._axes.text(center_x, center_y, "%s" % halo.id,
                    fontsize=self.font_size)

class HopParticleCallback(PlotCallback):
    """
    Adds particle positions for the members of each halo as identified
    by HOP. Along *axis* up to *max_number* groups in *hop_output* that are
    larger than *min_size* are plotted with *p_size* pixels per particle; 
    *alpha* determines the opacity of each particle.
    """
    def __init__(self, hop_output, axis, p_size=1.0,
                max_number=None, min_size=20, alpha=0.2):
        self.axis = axis
        self.hop_output = hop_output
        self.p_size = p_size
        self.max_number = max_number
        self.min_size = min_size
        self.alpha = alpha
    
    def __call__(self,plot):
        if self.max_number < 1: return
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        xf = lagos.axis_names[lagos.x_dict[plot.data.axis]]
        yf = lagos.axis_names[lagos.y_dict[plot.data.axis]]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        # now we loop over the haloes
        for halo in self.hop_output[:self.max_number]:
            size = halo.get_size()
            if size < self.min_size: continue
            colors = na.ones(size)
            plot._axes.hold(True)
            plot._axes.scatter(halo.get_positions(xf)*dx,
                halo.get_positions(yf)*dx, edgecolors="None",
                s=self.p_size, c='black', alpha=self.alpha)
            plot._axes.set_xlim(xx0,xx1)
            plot._axes.set_ylim(yy0,yy1)
            plot._axes.hold(False)

class FloorToValueInPlot(PlotCallback):
    def __init__(self):
        pass

    def __call__(self, plot):
        aa = plot.image._A
        min_val = aa[aa>0].min()
        aa[aa==0] = min_val


class VobozCircleCallback(PlotCallback):
    def __init__(self, voboz_output, axis, max_number=None,
                 annotate=False, min_size=20, font_size=8, print_halo_size=False):
        self.axis = axis
        self.voboz_output = voboz_output
        self.max_number = max_number
        self.annotate = annotate
        self.min_size = min_size
        self.font_size = font_size
        self.print_halo_size = print_halo_size

    def __call__(self, plot):
        from matplotlib.patches import Circle
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = _get_bounds(plot._axes.bbox)
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        for i,halo in enumerate(self.voboz_output[:self.max_number]):
            if (len(halo.particles) >= self.min_size):
                radius = halo.maximum_radius * dx
                center = halo.center_of_mass
                center_x = (center[xi] - x0)*dx
                center_y = (center[yi] - y0)*dy
                #print "voboz center = (%f,%f)" % (center[xi],center[yi])
                #print "voboz radius = %f" % halo.maximum_radius
                cir = Circle((center_x, center_y), radius, fill=False)
                plot._axes.add_patch(cir)
                if self.annotate:
                    if self.print_halo_size:
                        plot._axes.text(center_x, center_y, "%s" % len(halo.particles),
                        fontsize=self.font_size)
                    else:
                        plot._axes.text(center_x, center_y, "%s" % i,
                        fontsize=self.font_size)

class CoordAxesCallback(PlotCallback):
    """Creates x and y axes for a VMPlot. In the future, it will
    attempt to guess the proper units to use.

    """
    def __init__(self,unit=None,coords=False):
        PlotCallback.__init__(self)
        self.unit = unit
        self.coords = coords

    def __call__(self,plot):
        # 1. find out what the domain is
        # 2. pick a unit for it
        # 3. run self._axes.set_xlabel & self._axes.set_ylabel to actually lay shit down.
        # 4. adjust extent information to make sure labels are visable.

        # put the plot into data coordinates
        nx,ny = plot.image._A.shape
        dx = (plot.xlim[1] - plot.xlim[0])/nx
        dy = (plot.ylim[1] - plot.ylim[0])/ny

        unit_conversion = plot.data.hierarchy[plot.im["Unit"]]
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
        
        xlabel = '%s (%s)' % (lagos.axis_labels[plot.data.axis][0],plot.im["Unit"])
        ylabel = '%s (%s)' % (lagos.axis_labels[plot.data.axis][1],plot.im["Unit"])
        xticksize = nx/4.
        yticksize = ny/4.
        plot._axes.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([i*xticksize for i in range(0,5)]))
        plot._axes.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([i*yticksize for i in range(0,5)]))
        
        plot._axes.set_xlabel(xlabel,visible=True)
        plot._axes.set_ylabel(ylabel,visible=True)
        plot._figure.subplots_adjust(left=0.1,right=0.8)


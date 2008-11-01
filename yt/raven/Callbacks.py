"""
Callbacks to add additional functionality on to plots.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

import _MPL

class PlotCallback(object):
    def __init__(self, *args, **kwargs):
        pass

    def convert_to_pixels(self, plot, coord):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = plot._axes.bbox.get_bounds()
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        return ((coord[0] - x0)*dx, (coord[1] - y0)*dy)

class QuiverCallback(PlotCallback):
    def __init__(self, field_x, field_y, factor):
        """
        Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
        from the associated data, skipping every *factor* datapoints.
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
        pixX = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_x],
                             int(nx), int(ny),
                           (x0, x1, y0, y1),).transpose()
        pixY = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_y],
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
        xlim = na.logical_and(plot.data["px"] >= x0*0.9,
                              plot.data["px"] <= x1*1.1)
        ylim = na.logical_and(plot.data["py"] >= y0*0.9,
                              plot.data["py"] <= y1*1.1)
        wI = na.where(na.logical_and(xlim,ylim))
        xi, yi = na.mgrid[0:numPoints_x:numPoints_x/(self.factor*1j),\
                          0:numPoints_y:numPoints_y/(self.factor*1j)]
        x = (plot.data["px"][wI]-x0)*dx
        y = (plot.data["py"][wI]-y0)*dy
        z = plot.data[self.field][wI]
        if self.take_log: z=na.log10(z)
        zi = self.de.Triangulation(x,y).nn_interpolator(z)(xi,yi)
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
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        GLE = plot.data.gridLeftEdge
        GRE = plot.data.gridRightEdge
        px_index = lagos.x_dict[plot.data.axis]
        py_index = lagos.y_dict[plot.data.axis]
        left_edge_px = na.maximum((GLE[:,px_index]-x0)*dx, xx0)
        left_edge_py = na.maximum((GLE[:,py_index]-y0)*dy, yy0)
        right_edge_px = na.minimum((GRE[:,px_index]-x0)*dx, xx1)
        right_edge_py = na.minimum((GRE[:,py_index]-y0)*dy, yy1)
        print left_edge_px.min(), left_edge_px.max(), \
              right_edge_px.min(), right_edge_px.max(), \
              x0, x1, y0, y1
        verts = na.array(
                [(left_edge_px, left_edge_px, right_edge_px, right_edge_px),
                 (left_edge_py, right_edge_py, right_edge_py, left_edge_py)])
        visible =  ( right_edge_px - left_edge_px > self.min_pix ) & \
                   ( right_edge_px - left_edge_px > self.min_pix )
        verts=verts.transpose()[visible,:,:]
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
        l, b, width, height = plot._axes.bbox.get_bounds()
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
    def __init__(self, clumps, axis, plot_args = None):
        """
        Take a list of *clumps* and plot them as a set of contours.
        """
        self.clumps = clumps
        self.xf = lagos.axis_names[lagos.x_dict[axis]]
        self.yf = lagos.axis_names[lagos.y_dict[axis]]
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        
        nx, ny = plot.image._A.shape
        buff = na.zeros((nx,ny),dtype='float64')
        for i,clump in enumerate(reversed(self.clumps)):
            mylog.debug("Pixelizing contour %s", i)
            temp = _MPL.Pixelize(clump[self.xf],
                                 clump[self.yf],
                                 clump['dx'],
                                 clump['dx'],
                                 clump['dx']*0.0+i+1, # inits inside Pixelize
                                 int(nx), int(ny),
                             (x0, x1, y0, y1), 0).transpose()
            buff = na.maximum(temp, buff)
        self.rv = plot._axes.contour(buff, len(self.clumps)+1,
                                     **self.plot_args)
        plot._axes.hold(False)

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
        l, b, width, height = plot._axes.bbox.get_bounds()
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        radius = self.radius * dx
        center_x = (self.center[xi] - x0)*dx
        center_y = (self.center[yi] - y0)*dy
        cir = Circle((center_x, center_y), radius, fill=False,
                     **self.circle_args)
        plot._axes.add_patch(cir)
        if self.text is not None:
            plot._axes.text(center_x, center_y, "%s" % halo.id,
                            **self.text_args)

        

class HopCircleCallback(PlotCallback):
    def __init__(self, hop_output, axis, max_number=None,
                 annotate=False, min_size=20, font_size=8, print_halo_size=False):
        self.axis = axis
        self.hop_output = hop_output
        self.max_number = max_number
        self.annotate = annotate
        self.min_size = min_size
        self.font_size = font_size
        self.print_halo_size = print_halo_size

    def __call__(self, plot):
        from matplotlib.patches import Circle
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        l, b, width, height = plot._axes.bbox.get_bounds()
        xi = lagos.x_dict[plot.data.axis]
        yi = lagos.y_dict[plot.data.axis]
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        for halo in self.hop_output[:self.max_number]:
            radius = halo.maximum_radius() * dx
            center = halo.center_of_mass()
            center_x = (center[xi] - x0)*dx
            center_y = (center[yi] - y0)*dy
            cir = Circle((center_x, center_y), radius, fill=False)
            plot._axes.add_patch(cir)
            if self.annotate:
                if self.print_halo_size:
                    plot._axes.text(center_x, center_y, "%s" % len(halo.indices),
                    fontsize=self.font_size)
                else:
                    plot._axes.text(center_x, center_y, "%s" % halo.id,
                    fontsize=self.font_size)

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
        l, b, width, height = plot._axes.bbox.get_bounds()
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


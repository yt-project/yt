"""
Callbacks to add additional functionality on to plots.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
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

class QuiverCallback(PlotCallback):
    def __init__(self, field_x, field_y, factor):
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
        numPoints_x = plot.image._A.shape[0] / self.factor
        numPoints_y = plot.image._A.shape[1] / self.factor
        pixX = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_x],
                             int(numPoints_x), int(numPoints_y),
                           (x0, x1, y0, y1),).transpose()
        pixY = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             plot.data[self.field_y],
                             int(numPoints_x), int(numPoints_y),
                           (x0, x1, y0, y1),).transpose()
        X = na.mgrid[0:plot.image._A.shape[0]-1:numPoints_x*1j]# + 0.5*factor
        Y = na.mgrid[0:plot.image._A.shape[1]-1:numPoints_y*1j]# + 0.5*factor
        plot._axes.quiver(X,Y, pixX, -pixY)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ParticleCallback(PlotCallback):
    def __init__(self, axis, width, p_size=1.0, col='k'):
        PlotCallback.__init__(self)
        self.axis = axis
        self.width = width
        self.p_size = p_size
        self.color = col
        self.field_x = "particle_position_%s" % lagos.axis_names[lagos.x_dict[axis]]
        self.field_y = "particle_position_%s" % lagos.axis_names[lagos.y_dict[axis]]
        self.field_z = "particle_position_%s" % lagos.axis_names[axis]

    def __call__(self, plot):
        z0 = plot.data.center[self.axis] - self.width/2.0
        z1 = plot.data.center[self.axis] + self.width/2.0
        grids = plot.data._grids
        particles_x = na.concatenate([g[self.field_x] for g in grids]).ravel()
        particles_y = na.concatenate([g[self.field_y] for g in grids]).ravel()
        particles_z = na.concatenate([g[self.field_z] for g in grids]).ravel()
        if len(particles_x) == 0: return
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        # Now we rescale because our axes limits != data limits
        goodI = na.where( (particles_x < x1) & (particles_x > x0)
                        & (particles_y < y1) & (particles_y > y0)
                        & (particles_z < z1) & (particles_z > z0))
        particles_x = (particles_x[goodI] - x0) * (xx1-xx0)/(x1-x0)
        particles_y = (particles_y[goodI] - y0) * (yy1-yy0)/(y1-y0)
        plot._axes.hold(True)
        plot._axes.scatter(particles_x, particles_y, edgecolors='None',
                          s=self.p_size, c=self.color)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ContourCallback(PlotCallback):
    def __init__(self, field, ncont=5, factor=4, take_log=False, clim=None):
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
        plot._axes.contour(xi,yi,zi,self.ncont,colors='k')
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class GridBoundaryCallback(PlotCallback):
    def __init__(self, alpha=1.0, min_pix = 1):
        PlotCallback.__init__(self)
        self.alpha = alpha
        self.min_pix = min_pix

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        dx = plot.image._A.shape[0] / (x1-x0)
        dy = plot.image._A.shape[1] / (y1-y0)
        GLE = plot.data.gridLeftEdge
        GRE = plot.data.gridRightEdge
        px_index = lagos.x_dict[plot.data.axis]
        py_index = lagos.y_dict[plot.data.axis]
        left_edge_px = (GLE[:,px_index]-x0)*dx
        left_edge_py = (GLE[:,py_index]-y0)*dy
        right_edge_px = (GRE[:,px_index]-x0)*dx
        right_edge_py = (GRE[:,py_index]-y0)*dy
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
        min_exp_x = na.ceil(na.log10(w_min_x*plot.data.pf[unit])
                           /na.log10(self.factor))
        max_exp_x = na.floor(na.log10(w_max_x*plot.data.pf[unit])
                            /na.log10(self.factor))
        n_x = max_exp_x - min_exp_x + 1
        widths = na.logspace(min_exp_x, max_exp_x, num = n_x, base=self.factor)
        widths /= plot.data.pf[unit]
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
        PlotCallback.__init__(self)
        self.x = x
        self.y = y
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        plot._axes.hold(True)
        plot._axes.plot(self.x, self.y, **self.plot_args)
        plot._axes.hold(False)

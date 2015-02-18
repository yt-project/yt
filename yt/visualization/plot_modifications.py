"""

Callbacks to add additional functionality on to plots.



"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import h5py

from distutils.version import LooseVersion

from matplotlib.patches import Circle
from matplotlib.colors import colorConverter

from yt.funcs import *
from yt.extern.six import add_metaclass
from ._mpl_imports import *
from yt.utilities.physical_constants import \
    sec_per_Gyr, sec_per_Myr, \
    sec_per_kyr, sec_per_year, \
    sec_per_day, sec_per_hr
from yt.units.yt_array import YTQuantity, YTArray
from yt.visualization.image_writer import apply_colormap
from yt.utilities.lib.geometry_utils import triangle_plane_intersect

from . import _MPL

callback_registry = {}

class RegisteredCallback(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        callback_registry[name] = cls

@add_metaclass(RegisteredCallback)
class PlotCallback(object):
    def __init__(self, *args, **kwargs):
        pass

    def convert_to_plot(self, plot, coord, offset = True):
        # coord should be a 2 x ncoord array-like datatype.
        try:
            ncoord = np.array(coord).shape[1]
        except IndexError:
            ncoord = 1

        # Convert the data and plot limits to tiled numpy arrays so that
        # convert_to_plot is automatically vectorized.

        x0 = np.array(np.tile(plot.xlim[0],ncoord))
        x1 = np.array(np.tile(plot.xlim[1],ncoord))
        xx0 = np.tile(plot._axes.get_xlim()[0],ncoord)
        xx1 = np.tile(plot._axes.get_xlim()[1],ncoord)

        y0 = np.array(np.tile(plot.ylim[0],ncoord))
        y1 = np.array(np.tile(plot.ylim[1],ncoord))
        yy0 = np.tile(plot._axes.get_ylim()[0],ncoord)
        yy1 = np.tile(plot._axes.get_ylim()[1],ncoord)

        ccoord = np.array(coord)

        # We need a special case for when we are only given one coordinate.
        if ccoord.shape == (2,):
            return ((ccoord[0]-x0)/(x1-x0)*(xx1-xx0) + xx0,
                    (ccoord[1]-y0)/(y1-y0)*(yy1-yy0) + yy0)
        else:
            return ((ccoord[0][:]-x0)/(x1-x0)*(xx1-xx0) + xx0,
                    (ccoord[1][:]-y0)/(y1-y0)*(yy1-yy0) + yy0)

    def pixel_scale(self, plot):
        x0, x1 = np.array(plot.xlim)
        xx0, xx1 = plot._axes.get_xlim()
        dx = (xx1 - xx0)/(x1 - x0)

        y0, y1 = np.array(plot.ylim)
        yy0, yy1 = plot._axes.get_ylim()
        dy = (yy1 - yy0)/(y1 - y0)

        return (dx,dy)

    def _set_font_properties(self, plot, labels, **kwargs):
        """
        This sets all of the text instances created by a callback to have
        the same font size and properties as all of the other fonts in the
        figure.  If kwargs are set, they trump the defaults.
        """
        # This is a little messy because there is no trivial way to update
        # a MPL.font_manager.FontProperties object with new attributes
        # aside from setting them individually.  So we pick out the relevant
        # MPL.Text() kwargs from the local kwargs and let them trump the
        # defaults.
        local_font_properties = plot.font_properties.copy()

        # Turn off the default TT font file, otherwise none of this works.
        local_font_properties.set_file(None)
        local_font_properties.set_family('stixgeneral')

        if kwargs.has_key('family'): 
            local_font_properties.set_family(kwargs['family'])
        if kwargs.has_key('file'): 
            local_font_properties.set_file(kwargs['file'])
        if kwargs.has_key('fontconfig_pattern'): 
            local_font_properties.set_fontconfig_pattern(kwargs['fontconfig_pattern'])
        if kwargs.has_key('name'): 
            local_font_properties.set_name(kwargs['name'])
        if kwargs.has_key('size'): 
            local_font_properties.set_size(kwargs['size'])
        if kwargs.has_key('slant'): 
            local_font_properties.set_slant(kwargs['slant'])
        if kwargs.has_key('stretch'): 
            local_font_properties.set_stretch(kwargs['stretch'])
        if kwargs.has_key('style'): 
            local_font_properties.set_style(kwargs['style'])
        if kwargs.has_key('variant'): 
            local_font_properties.set_variant(kwargs['variant'])
        if kwargs.has_key('weight'): 
            local_font_properties.set_weight(kwargs['weight'])

        # For each label, set the font properties and color to the figure
        # defaults if not already set in the callback itself
        for label in labels:
            if plot.font_color is not None and not kwargs.has_key('color'):
                label.set_color(plot.font_color)
            label.set_fontproperties(local_font_properties)

class VelocityCallback(PlotCallback):
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
    _type_name = "velocity"
    def __init__(self, factor=16, scale=None, scale_units=None, normalize=False):
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale  = scale
        self.scale_units = scale_units
        self.normalize = normalize

    def __call__(self, plot):
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            qcb = CuttingQuiverCallback("cutting_plane_velocity_x",
                                        "cutting_plane_velocity_y",
                                        self.factor)
        else:
            ax = plot.data.axis
            (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                        plot.data.ds.coordinates.y_axis[ax])
            axis_names = plot.data.ds.coordinates.axis_name
            xv = "velocity_%s" % (axis_names[xi])
            yv = "velocity_%s" % (axis_names[yi])

            bv = plot.data.get_field_parameter("bulk_velocity")
            if bv is not None:
                bv_x = bv[xi]
                bv_y = bv[yi]
            else: bv_x = bv_y = YTQuantity(0, 'cm/s')

            qcb = QuiverCallback(xv, yv, self.factor, scale=self.scale, 
                                 scale_units=self.scale_units, 
                                 normalize=self.normalize, bv_x=bv_x, bv_y=bv_y)
        return qcb(plot)

class MagFieldCallback(PlotCallback):
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
    _type_name = "magnetic_field"
    def __init__(self, factor=16, scale=None, scale_units=None, normalize=False):
        PlotCallback.__init__(self)
        self.factor = factor
        self.scale  = scale
        self.scale_units = scale_units
        self.normalize = normalize

    def __call__(self, plot):
        # Instantiation of these is cheap
        if plot._type_name == "CuttingPlane":
            qcb = CuttingQuiverCallback("cutting_plane_bx",
                                        "cutting_plane_by",
                                        self.factor)
        else:
            xax = plot.data.ds.coordinates.x_axis[plot.data.axis]
            yax = plot.data.ds.coordinates.y_axis[plot.data.axis]
            axis_names = plot.data.ds.coordinates.axis_name
            xv = "magnetic_field_%s" % (axis_names[xax])
            yv = "magnetic_field_%s" % (axis_names[yax])
            qcb = QuiverCallback(xv, yv, self.factor, scale=self.scale, scale_units=self.scale_units, normalize=self.normalize)
        return qcb(plot)

class QuiverCallback(PlotCallback):
    """
    annotate_quiver(field_x, field_y, factor=16, scale=None, scale_units=None, 
                    normalize=False, bv_x=0, bv_y=0):

    Adds a 'quiver' plot to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints
    *scale* is the data units per arrow length unit using *scale_units* 
    (see matplotlib.axes.Axes.quiver for more info)
    """
    _type_name = "quiver"
    def __init__(self, field_x, field_y, factor=16, scale=None, scale_units=None, normalize=False, bv_x=0, bv_y=0):
        PlotCallback.__init__(self)
        self.field_x = field_x
        self.field_y = field_y
        self.bv_x = bv_x
        self.bv_y = bv_y
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
        # periodicity
        ax = plot.data.axis
        ds = plot.data.ds
        (xi, yi) = (ds.coordinates.x_axis[ax],
                    ds.coordinates.y_axis[ax])
        period_x = ds.domain_width[xi]
        period_y = ds.domain_width[yi]
        periodic = int(any(ds.periodicity))
        fv_x = plot.data[self.field_x]
        if self.bv_x != 0.0:
            # Workaround for 0.0 without units
            fv_x -= self.bv_x
        fv_y = plot.data[self.field_y]
        if self.bv_y != 0.0:
            # Workaround for 0.0 without units
            fv_y -= self.bv_y
        pixX = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             fv_x,
                             int(nx), int(ny),
                             (x0, x1, y0, y1), 0, # bounds, antialias
                             (period_x, period_y), periodic,
                           ).transpose()
        pixY = _MPL.Pixelize(plot.data['px'],
                             plot.data['py'],
                             plot.data['pdx'],
                             plot.data['pdy'],
                             fv_y,
                             int(nx), int(ny),
                             (x0, x1, y0, y1), 0, # bounds, antialias
                             (period_x, period_y), periodic,
                           ).transpose()
        X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                          np.linspace(yy0,yy1,ny,endpoint=True))
        if self.normalize:
            nn = np.sqrt(pixX**2 + pixY**2)
            pixX /= nn
            pixY /= nn
        plot._axes.quiver(X,Y, pixX, pixY, scale=self.scale, scale_units=self.scale_units)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ContourCallback(PlotCallback):
    """
    annotate_contour(field, ncont=5, factor=4, take_log=None, clim=None,
                     plot_args=None, label=False, label_args=None,
                     data_source=None):

    Add contours in *field* to the plot.  *ncont* governs the number of
    contours generated, *factor* governs the number of points used in the
    interpolation, *take_log* governs how it is contoured and *clim* gives
    the (upper, lower) limits for contouring.  An alternate data source can be
    specified with *data_source*, but by default the plot's data source will be
    queried.
    """
    _type_name = "contour"
    def __init__(self, field, ncont=5, factor=4, clim=None,
                 plot_args = None, label = False, take_log = None, 
                 label_args = None, data_source = None):
        PlotCallback.__init__(self)
        self.ncont = ncont
        self.field = field
        self.factor = factor
        self.clim = clim
        self.take_log = take_log
        if plot_args is None: plot_args = {'colors':'k'}
        self.plot_args = plot_args
        self.label = label
        if label_args is None:
            label_args = {}
        self.label_args = label_args
        self.data_source = data_source

    def __call__(self, plot):
        # These need to be in code_length
        x0, x1 = (v.in_units("code_length") for v in plot.xlim)
        y0, y1 = (v.in_units("code_length") for v in plot.ylim)

        # These are in plot coordinates, which may not be code coordinates.
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()

        plot._axes.hold(True)
        
        numPoints_x = plot.image._A.shape[0]
        numPoints_y = plot.image._A.shape[1]
        
        # Multiply by dx and dy to go from data->plot
        dx = (xx1 - xx0) / (x1-x0)
        dy = (yy1 - yy0) / (y1-y0)

        # We want xi, yi in plot coordinates
        xi, yi = np.mgrid[xx0:xx1:numPoints_x/(self.factor*1j),
                          yy0:yy1:numPoints_y/(self.factor*1j)]
        data = self.data_source or plot.data

        if plot._type_name in ['CuttingPlane','Projection','Slice']:
            if plot._type_name == 'CuttingPlane':
                x = data["px"]*dx
                y = data["py"]*dy
                z = data[self.field]
            elif plot._type_name in ['Projection','Slice']:
                #Makes a copy of the position fields "px" and "py" and adds the
                #appropriate shift to the copied field.  

                AllX = np.zeros(data["px"].size, dtype='bool')
                AllY = np.zeros(data["py"].size, dtype='bool')
                XShifted = data["px"].copy()
                YShifted = data["py"].copy()
                dom_x, dom_y = plot._period
                for shift in np.mgrid[-1:1:3j]:
                    xlim = ((data["px"] + shift*dom_x >= x0) &
                            (data["px"] + shift*dom_x <= x1))
                    ylim = ((data["py"] + shift*dom_y >= y0) &
                            (data["py"] + shift*dom_y <= y1))
                    XShifted[xlim] += shift * dom_x
                    YShifted[ylim] += shift * dom_y
                    AllX |= xlim
                    AllY |= ylim
            
                # At this point XShifted and YShifted are the shifted arrays of
                # position data in data coordinates
                wI = (AllX & AllY)

                # This converts XShifted and YShifted into plot coordinates
                x = ((XShifted[wI]-x0)*dx).ndarray_view() + xx0
                y = ((YShifted[wI]-y0)*dy).ndarray_view() + yy0
                z = data[self.field][wI]
        
            # Both the input and output from the triangulator are in plot
            # coordinates
            if LooseVersion(matplotlib.__version__) < LooseVersion("1.4.0"):
                from matplotlib.delaunay.triangulate import Triangulation as \
                    triang
                zi = triang(x,y).nn_interpolator(z)(xi,yi)
            else:
                from matplotlib.tri import Triangulation, LinearTriInterpolator
                triangulation = Triangulation(x, y)
                zi = LinearTriInterpolator(triangulation, z)(xi,yi)
        elif plot._type_name == 'OffAxisProjection':
            zi = plot.frb[self.field][::self.factor,::self.factor].transpose()
        
        if self.take_log is None:
            field = data._determine_fields([self.field])[0]
            self.take_log = plot.ds._get_field_info(*field).take_log

        if self.take_log: zi=np.log10(zi)

        if self.take_log and self.clim is not None: 
            self.clim = (np.log10(self.clim[0]), np.log10(self.clim[1]))
        
        if self.clim is not None: 
            self.ncont = np.linspace(self.clim[0], self.clim[1], self.ncont)
        
        cset = plot._axes.contour(xi,yi,zi,self.ncont, **self.plot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)
        
        if self.label:
            plot._axes.clabel(cset, **self.label_args)
        

class GridBoundaryCallback(PlotCallback):
    """
    annotate_grids(alpha=0.7, min_pix=1, min_pix_ids=20, draw_ids=False,
                   periodic=True, min_level=None, max_level=None,
                   cmap='B-W LINEAR_r', edgecolors=None, linewidth=1.0):

    Draws grids on an existing PlotWindow object.  Adds grid boundaries to a
    plot, optionally with alpha-blending. By default, colors different levels of
    grids with different colors going from white to black, but you can change to
    any arbitrary colormap with cmap keyword, to all black grid edges for all
    levels with cmap=None and edgecolors=None, or to an arbitrary single color
    for grid edges with edgecolors='YourChosenColor' defined in any of the
    standard ways (e.g., edgecolors='white', edgecolors='r',
    edgecolors='#00FFFF', or edgecolor='0.3', where the last is a float in 0-1
    scale indicating gray).  Note that setting edgecolors overrides cmap if you
    have both set to non-None values.  Cutoff for display is at min_pix
    wide. draw_ids puts the grid id in the corner of the grid.  (Not so great in
    projections...).  One can set min and maximum level of grids to display, and
    can change the linewidth of the displayed grids.
    """
    _type_name = "grids"

    def __init__(self, alpha=0.7, min_pix=1, min_pix_ids=20, draw_ids=False,
                 periodic=True, min_level=None, max_level=None,
                 cmap='B-W LINEAR_r', edgecolors=None, linewidth=1.0):
        PlotCallback.__init__(self)
        self.alpha = alpha
        self.min_pix = min_pix
        self.min_pix_ids = min_pix_ids
        self.draw_ids = draw_ids  # put grid numbers in the corner.
        self.periodic = periodic
        self.min_level = min_level
        self.max_level = max_level
        self.linewidth = linewidth
        self.cmap = cmap
        self.edgecolors = edgecolors

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        (dx, dy) = self.pixel_scale(plot)
        (xpix, ypix) = plot.image._A.shape
        ax = plot.data.axis
        px_index = plot.data.ds.coordinates.x_axis[ax]
        py_index = plot.data.ds.coordinates.y_axis[ax]
        DW = plot.data.ds.domain_width
        if self.periodic:
            pxs, pys = np.mgrid[-1:1:3j,-1:1:3j]
        else:
            pxs, pys = np.mgrid[0:0:1j,0:0:1j]
        GLE, GRE, levels = [], [], []
        for block, mask in plot.data.blocks:
            GLE.append(block.LeftEdge.in_units("code_length"))
            GRE.append(block.RightEdge.in_units("code_length"))
            levels.append(block.Level)
        if len(GLE) == 0: return
        # Retain both units and registry
        GLE = YTArray(GLE, input_units = GLE[0].units)
        GRE = YTArray(GRE, input_units = GRE[0].units)
        levels = np.array(levels)
        min_level = self.min_level or 0
        max_level = self.max_level or levels.max()

        # sorts the three arrays in order of ascending level - this makes images look nicer
        new_indices = np.argsort(levels)
        levels = levels[new_indices]
        GLE = GLE[new_indices]
        GRE = GRE[new_indices]
        
        for px_off, py_off in zip(pxs.ravel(), pys.ravel()):
            pxo = px_off * DW[px_index]
            pyo = py_off * DW[py_index]
            left_edge_x = np.array((GLE[:,px_index]+pxo-x0)*dx) + xx0
            left_edge_y = np.array((GLE[:,py_index]+pyo-y0)*dy) + yy0
            right_edge_x = np.array((GRE[:,px_index]+pxo-x0)*dx) + xx0
            right_edge_y = np.array((GRE[:,py_index]+pyo-y0)*dy) + yy0
            xwidth = xpix * (right_edge_x - left_edge_x) / (xx1 - xx0)
            ywidth = ypix * (right_edge_y - left_edge_y) / (yy1 - yy0)
            visible = np.logical_and(
                np.logical_and(xwidth > self.min_pix, ywidth > self.min_pix),
                np.logical_and(levels >= min_level, levels <= max_level))

            # Grids can either be set by edgecolors OR a colormap.
            if self.edgecolors is not None:
                edgecolors = colorConverter.to_rgba(
                    self.edgecolors, alpha=self.alpha)
            else:  # use colormap if not explicity overridden by edgecolors
                if self.cmap is not None:
                    color_bounds = [0,plot.data.ds.index.max_level]
                    edgecolors = apply_colormap(
                        levels[visible]*1.0, color_bounds=color_bounds,
                        cmap_name=self.cmap)[0,:,:]*1.0/255.
                    edgecolors[:,3] = self.alpha
                else:
                    edgecolors = (0.0,0.0,0.0,self.alpha)

            if visible.nonzero()[0].size == 0: continue
            verts = np.array(
                [(left_edge_x, left_edge_x, right_edge_x, right_edge_x),
                 (left_edge_y, right_edge_y, right_edge_y, left_edge_y)])
            verts=verts.transpose()[visible,:,:]
            grid_collection = matplotlib.collections.PolyCollection(
                verts, facecolors="none", edgecolors=edgecolors,
                linewidth=self.linewidth)
            plot._axes.hold(True)
            plot._axes.add_collection(grid_collection)

            if self.draw_ids:
                visible_ids = np.logical_and(
                    np.logical_and(xwidth > self.min_pix_ids,
                                   ywidth > self.min_pix_ids),
                    np.logical_and(levels >= min_level, levels <= max_level))
                active_ids = np.unique(plot.data['grid_indices'])
                for i in np.where(visible_ids)[0]:
                    plot._axes.text(
                        left_edge_x[i] + (2 * (xx1 - xx0) / xpix),
                        left_edge_y[i] + (2 * (yy1 - yy0) / ypix),
                        "%d" % active_ids[i], clip_on=True)
            plot._axes.hold(False)

class StreamlineCallback(PlotCallback):
    """
    annotate_streamlines(field_x, field_y, factor = 16,
                         density = 1, plot_args=None):

    Add streamlines to any plot, using the *field_x* and *field_y*
    from the associated data, skipping every *factor* datapoints like
    'quiver'. *density* is the index of the amount of the streamlines.
    """
    _type_name = "streamlines"
    def __init__(self, field_x, field_y, factor = 16,
                 density = 1, plot_args=None):
        PlotCallback.__init__(self)
        self.field_x = field_x
        self.field_y = field_y
        self.factor = factor
        self.dens = density
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args
        
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
        X,Y = (np.linspace(xx0,xx1,nx,endpoint=True),
               np.linspace(yy0,yy1,ny,endpoint=True))
        streamplot_args = {'x': X, 'y': Y, 'u':pixX, 'v': pixY,
                           'density': self.dens}
        streamplot_args.update(self.plot_args)
        plot._axes.streamplot(**streamplot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class LabelCallback(PlotCallback):
    """
    This adds a label to the plot.
    """
    _type_name = "axis_label"
    def __init__(self, label):
        PlotCallback.__init__(self)
        self.label = label

    def __call__(self, plot):
        plot._figure.subplots_adjust(hspace=0, wspace=0,
                                     bottom=0.1, top=0.9,
                                     left=0.0, right=1.0)
        plot._axes.set_xlabel(self.label)
        plot._axes.set_ylabel(self.label)

class LinePlotCallback(PlotCallback):
    """
    annotate_line(x, y, plot_args = None)

    Over plot *x* and *y* with *plot_args* fed into the plot.
    """
    _type_name = "line"
    def __init__(self, x, y, plot_args = None):
        PlotCallback.__init__(self)
        self.x = x
        self.y = y
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        plot._axes.hold(True)
        plot._axes.plot(self.x, self.y, **self.plot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ImageLineCallback(LinePlotCallback):
    """
    annotate_image_line(p1, p2, data_coords=False, plot_args = None)

    Plot from *p1* to *p2* (image plane coordinates)
    with *plot_args* fed into the plot.
    """
    _type_name = "image_line"
    _plot_args = {'color':'white', 'linewidth':2}
    def __init__(self, p1, p2, data_coords=False, plot_args = None):
        PlotCallback.__init__(self)
        self.p1 = p1
        self.p2 = p2
        if plot_args is None: plot_args = self._plot_args
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
    """
    annotate_cquiver(field_x, field_y, factor)

    Get a quiver plot on top of a cutting plane, using *field_x* and
    *field_y*, skipping every *factor* datapoint in the discretization.
    """
    _type_name = "cquiver"
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
        nx = plot.image._A.shape[0] / self.factor
        ny = plot.image._A.shape[1] / self.factor
        indices = np.argsort(plot.data['dx'])[::-1]
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
        X,Y = np.meshgrid(np.linspace(xx0,xx1,nx,endpoint=True),
                          np.linspace(yy0,yy1,ny,endpoint=True))
        plot._axes.quiver(X,Y, pixX, pixY)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class ClumpContourCallback(PlotCallback):
    """
    annotate_clumps(clumps, plot_args = None)

    Take a list of *clumps* and plot them as a set of contours.
    """
    _type_name = "clumps"
    def __init__(self, clumps, plot_args = None):
        self.clumps = clumps
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()

        extent = [xx0,xx1,yy0,yy1]

        plot._axes.hold(True)

        ax = plot.data.axis
        px_index = plot.data.ds.coordinates.x_axis[ax]
        py_index = plot.data.ds.coordinates.y_axis[ax]

        xf = plot.data.ds.coordinates.axis_name[px_index]
        yf = plot.data.ds.coordinates.axis_name[py_index]
        dxf = "d%s" % xf
        dyf = "d%s" % yf

        nx, ny = plot.image._A.shape
        buff = np.zeros((nx,ny),dtype='float64')
        for i,clump in enumerate(reversed(self.clumps)):
            mylog.info("Pixelizing contour %s", i)

            xf_copy = clump[xf].copy().in_units("code_length")
            yf_copy = clump[yf].copy().in_units("code_length")

            temp = _MPL.Pixelize(xf_copy, yf_copy,
                                 clump[dxf].in_units("code_length")/2.0,
                                 clump[dyf].in_units("code_length")/2.0,
                                 clump[dxf].d*0.0+i+1, # inits inside Pixelize
                                 int(nx), int(ny),
                             (x0, x1, y0, y1), 0).transpose()
            buff = np.maximum(temp, buff)
        self.rv = plot._axes.contour(buff, np.unique(buff),
                                     extent=extent, **self.plot_args)
        plot._axes.hold(False)

class ArrowCallback(PlotCallback):
    """
    annotate_arrow(pos, code_size, plot_args = None)

    This adds an arrow pointing at *pos* with size *code_size* in code
    units.  *plot_args* is a dict fed to matplotlib with arrow properties.
    """
    _type_name = "arrow"
    def __init__(self, pos, code_size, plot_args = None):
        self.pos = pos
        if isinstance(code_size, YTArray):
            code_size = code_size.in_units('code_length')
        if not iterable(code_size):
            code_size = (code_size, code_size)
        self.code_size = code_size
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        if len(self.pos) == 3:
            ax = plot.data.axis
            (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                        plot.data.ds.coordinates.y_axis[ax])
            pos = self.pos[xi], self.pos[yi]
        else: pos = self.pos
        if isinstance(self.code_size[1], basestring):
            code_size = plot.data.ds.quan(*self.code_size)
            code_size = code_size.in_units('code_length').value
            self.code_size = (code_size, code_size)
        from matplotlib.patches import Arrow
        # Now convert the pixels to code information
        x, y = self.convert_to_plot(plot, pos)
        x1, y1 = pos[0]+self.code_size[0], pos[1]+self.code_size[1]
        x1, y1 = self.convert_to_plot(plot, (x1, y1), False)
        dx, dy = x1 - x, y1 - y
        arrow = Arrow(x-dx, y-dy, dx, dy, **self.plot_args)
        plot._axes.add_patch(arrow)

class PointAnnotateCallback(PlotCallback):
    """
    annotate_point(pos, text, text_args = None)

    This adds *text* at position *pos*, where *pos* is in code-space.
    *text_args* is a dict fed to the text placement code.
    """
    _type_name = "point"
    def __init__(self, pos, text, text_args = None):
        self.pos = pos
        self.text = text
        if text_args is None: text_args = {}
        self.text_args = text_args

    def __call__(self, plot):
        if len(self.pos) == 3:
            ax = plot.data.axis
            (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                        plot.data.ds.coordinates.y_axis[ax])
            pos = self.pos[xi], self.pos[yi]
        else: pos = self.pos
        x,y = self.convert_to_plot(plot, pos)
        
        # Set the font properties of text from this callback to be
        # consistent with other text labels in this figure
        label = plot._axes.text(x, y, self.text, **self.text_args)
        self._set_font_properties(plot, [label], **self.text_args)

class MarkerAnnotateCallback(PlotCallback):
    """
    annotate_marker(pos, marker='x', plot_args=None)

    Adds text *marker* at *pos* in code units.  *plot_args* is a dict
    that will be forwarded to the plot command.
    """
    _type_name = "marker"
    def __init__(self, pos, marker='x', plot_args=None):
        self.pos = pos
        self.marker = marker
        if plot_args is None: plot_args = {}
        self.plot_args = plot_args

    def __call__(self, plot):
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        if len(self.pos) == 3:
            ax = plot.data.axis
            (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                        plot.data.ds.coordinates.y_axis[ax])
            pos = self.pos[xi], self.pos[yi]
        elif len(self.pos) == 2:
            pos = self.pos
        x,y = self.convert_to_plot(plot, pos)
        plot._axes.hold(True)
        plot._axes.scatter(x,y, marker = self.marker, **self.plot_args)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

class SphereCallback(PlotCallback):
    """
    annotate_sphere(center, radius, circle_args = None,
                    text = None, text_args = None)
    
    A sphere centered at *center* in code units with radius *radius* in
    code units will be created, with optional *circle_args*, *text*, and
    *text_args*.
    """
    _type_name = "sphere"
    def __init__(self, center, radius, circle_args = None,
                 text = None, text_args = None):
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

        if iterable(self.radius):
            self.radius = plot.data.ds.quan(self.radius[0], self.radius[1])
            self.radius = np.float64(self.radius.in_units(plot.xlim[0].units))

        radius = self.radius * self.pixel_scale(plot)[0]

        if plot.data.axis == 4:
            (xi, yi) = (0, 1)
        else:
            ax = plot.data.axis
            (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                        plot.data.ds.coordinates.y_axis[ax])

        (center_x,center_y) = self.convert_to_plot(plot,(self.center[xi], self.center[yi]))
        
        cir = Circle((center_x, center_y), radius, **self.circle_args)
        plot._axes.add_patch(cir)
        if self.text is not None:
            plot._axes.text(center_x, center_y, self.text,
                            **self.text_args)


class TextLabelCallback(PlotCallback):
    """
    annotate_text(pos, text, data_coords=False, text_args=None, bbox_args=None)

    Accepts a position in (0..1, 0..1) of the image, some text and
    optionally some text arguments. If data_coords is True,
    position will be in code units instead of image coordinates.  If you desire
    a bounding box around your text, set one with the bbox_args dictionary 
    keyword.
    """
    _type_name = "text"
    def __init__(self, pos, text, data_coords=False, text_args=None, 
                 bbox_args=None):
        self.pos = pos
        self.text = text
        self.data_coords = data_coords
        if text_args is None: text_args = {}
        self.text_args = text_args
        if bbox_args is None: bbox_args = {}
        self.bbox_args = bbox_args

    def __call__(self, plot):
        kwargs = self.text_args.copy()
        if self.data_coords and len(plot.image._A.shape) == 2:
            if len(self.pos) == 3:
                ax = plot.data.axis
                (xi, yi) = (plot.data.ds.coordinates.x_axis[ax],
                            plot.data.ds.coordinates.y_axis[ax])
                pos = self.pos[xi], self.pos[yi]
            else: pos = self.pos
            x,y = self.convert_to_plot(plot, pos)
        else:
            x, y = self.pos
            if not self.data_coords:
                kwargs["transform"] = plot._axes.transAxes

        # Set the font properties of text from this callback to be
        # consistent with other text labels in this figure
        label = plot._axes.text(x, y, self.text, bbox=self.bbox_args, **kwargs)
        self._set_font_properties(plot, [label], **kwargs)

class HaloCatalogCallback(PlotCallback):
    """
    annotate_halos(halo_catalog, circle_kwargs=None,
        width=None, annotate_field=None,
        font_kwargs=None, factor = 1.0)

    Plots circles at the locations of all the halos
    in a halo catalog with radii corresponding to the
    virial radius of each halo. 

    circle_kwargs: Contains the arguments controlling the
        appearance of the circles, supplied to the 
        Matplotlib patch Circle.
    width: the width over which to select halos to plot,
        useful when overplotting to a slice plot. Accepts
        a tuple in the form (1.0, 'Mpc').
    annotate_field: Accepts a field contained in the 
        halo catalog to add text to the plot near the halo.
        Example: annotate_field = 'particle_mass' will
        write the halo mass next to each halo.
    font_kwargs: Contains the arguments controlling the text
        appearance of the annotated field.
    factor: A number the virial radius is multiplied by for
        plotting the circles. Ex: factor = 2.0 will plot
        circles with twice the radius of each halo virial radius.
    """

    _type_name = 'halos'
    region = None
    _descriptor = None

    def __init__(self, halo_catalog, circle_kwargs=None, 
            width=None, annotate_field=None,
            font_kwargs=None, factor = 1.0):

        PlotCallback.__init__(self)
        self.halo_catalog = halo_catalog
        self.width = width
        self.annotate_field = annotate_field
        if font_kwargs is None:
            font_kwargs = {'color':'white'}
        self.font_kwargs = font_kwargs
        self.factor = factor
        if circle_kwargs is None:
            circle_kwargs = {'edgecolor':'white', 'facecolor':'None'}
        self.circle_kwargs = circle_kwargs

    def __call__(self, plot):
        data = plot.data
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        
        halo_data= self.halo_catalog.halos_ds.all_data()
        axis_names = plot.data.ds.coordinates.axis_name
        xax = plot.data.ds.coordinates.x_axis[data.axis]
        yax = plot.data.ds.coordinates.y_axis[data.axis]
        field_x = "particle_position_%s" % axis_names[xax]
        field_y = "particle_position_%s" % axis_names[yax]
        field_z = "particle_position_%s" % axis_names[data.axis]
        plot._axes.hold(True)

        # Set up scales for pixel size and original data
        pixel_scale = self.pixel_scale(plot)[0]
        data_scale = data.ds.length_unit
        units = data_scale.units

        # Convert halo positions to code units of the plotted data
        # and then to units of the plotted window
        px = halo_data[field_x][:].in_units(units) / data_scale
        py = halo_data[field_y][:].in_units(units) / data_scale
        px, py = self.convert_to_plot(plot,[px,py])

        # Convert halo radii to a radius in pixels
        radius = halo_data['virial_radius'][:].in_units(units)
        radius = np.array(radius*pixel_scale*self.factor/data_scale)
        
        if self.width:
            pz = halo_data[field_z][:].in_units(units)/data_scale
            pz = data.ds.arr(pz, 'code_length')
            c = data.center[data.axis]

            # I should catch an error here if width isn't in this form
            # but I dont really want to reimplement get_sanitized_width...
            width = data.ds.arr(self.width[0], self.width[1]).in_units('code_length')

            indices = np.where((pz > c-width) & (pz < c+width))

            px = px[indices]
            py = py[indices]
            radius = radius[indices]

        for x,y,r in zip(px, py, radius):
            plot._axes.add_artist(Circle(xy=(x,y), 
                radius = r, **self.circle_kwargs)) 

        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)

        if self.annotate_field:
            annotate_dat = halo_data[self.annotate_field]
            texts = ['{:g}'.format(float(dat))for dat in annotate_dat]
            labels = []
            for pos_x, pos_y, t in zip(px, py, texts): 
                labels.append(plot._axes.text(pos_x, pos_y, t, **self.font_kwargs))

            # Set the font properties of text from this callback to be
            # consistent with other text labels in this figure
            self._set_font_properties(plot, labels, **self.font_kwargs)

class ParticleCallback(PlotCallback):
    """
    annotate_particles(width, p_size=1.0, col='k', marker='o', stride=1.0,
                       ptype=None, minimum_mass=None, alpha=1.0)

    Adds particle positions, based on a thick slab along *axis* with a
    *width* along the line of sight.  *p_size* controls the number of
    pixels per particle, and *col* governs the color.  *ptype* will
    restrict plotted particles to only those that are of a given type.
    Particles with masses below *minimum_mass* will not be plotted.
    *alpha* determines the opacity of the marker symbol used in the scatter
    plot.
    """
    _type_name = "particles"
    region = None
    _descriptor = None
    def __init__(self, width, p_size=1.0, col='k', marker='o', stride=1.0,
                 ptype='all', minimum_mass=None, alpha=1.0):
        PlotCallback.__init__(self)
        self.width = width
        self.p_size = p_size
        self.color = col
        self.marker = marker
        self.stride = stride
        self.ptype = ptype
        self.minimum_mass = minimum_mass
        self.alpha = alpha

    def __call__(self, plot):
        data = plot.data
        if iterable(self.width):
            self.width = np.float64(plot.data.ds.quan(self.width[0], self.width[1]))
        # we construct a recantangular prism
        x0, x1 = plot.xlim
        y0, y1 = plot.ylim
        xx0, xx1 = plot._axes.get_xlim()
        yy0, yy1 = plot._axes.get_ylim()
        reg = self._get_region((x0,x1), (y0,y1), plot.data.axis, data)
        ax = data.axis
        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]
        axis_names = plot.data.ds.coordinates.axis_name
        field_x = "particle_position_%s" % axis_names[xax]
        field_y = "particle_position_%s" % axis_names[yax]
        pt = self.ptype
        gg = ( ( reg[pt, field_x] >= x0 ) & ( reg[pt, field_x] <= x1 )
           &   ( reg[pt, field_y] >= y0 ) & ( reg[pt, field_y] <= y1 ) )
        if self.minimum_mass is not None:
            gg &= (reg[pt, "particle_mass"] >= self.minimum_mass)
            if gg.sum() == 0: return
        plot._axes.hold(True)
        px, py = self.convert_to_plot(plot,
                    [np.array(reg[pt, field_x][gg][::self.stride]),
                     np.array(reg[pt, field_y][gg][::self.stride])])
        plot._axes.scatter(px, py, edgecolors='None', marker=self.marker,
                           s=self.p_size, c=self.color,alpha=self.alpha)
        plot._axes.set_xlim(xx0,xx1)
        plot._axes.set_ylim(yy0,yy1)
        plot._axes.hold(False)


    def _get_region(self, xlim, ylim, axis, data):
        LE, RE = [None]*3, [None]*3
        ds = data.ds
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        zax = axis
        LE[xax], RE[xax] = xlim
        LE[yax], RE[yax] = ylim
        LE[zax] = data.center[zax].ndarray_view() - self.width*0.5
        RE[zax] = data.center[zax].ndarray_view() + self.width*0.5
        if self.region is not None \
            and np.all(self.region.left_edge <= LE) \
            and np.all(self.region.right_edge >= RE):
            return self.region
        self.region = data.ds.region(data.center, LE, RE)
        return self.region

class TitleCallback(PlotCallback):
    """
    annotate_title(title)

    Accepts a *title* and adds it to the plot
    """
    _type_name = "title"
    def __init__(self, title):
        PlotCallback.__init__(self)
        self.title = title

    def __call__(self,plot):
        plot._axes.set_title(self.title)
        # Set the font properties of text from this callback to be
        # consistent with other text labels in this figure
        label = plot._axes.title
        self._set_font_properties(plot, [label])

class TriangleFacetsCallback(PlotCallback):
    """ 
    annotate_triangle_facets(triangle_vertices, plot_args=None )

    Intended for representing a slice of a triangular faceted 
    geometry in a slice plot. 

    Uses a set of *triangle_vertices* to find all trangles the plane of a 
    SlicePlot intersects with. The lines between the intersection points 
    of the triangles are then added to the plot to create an outline
    of the geometry represented by the triangles. 
    """
    _type_name = "triangle_facets"
    def __init__(self, triangle_vertices, plot_args=None):
        super(TriangleFacetsCallback, self).__init__()
        self.plot_args = {} if plot_args is None else plot_args
        self.vertices = triangle_vertices

    def __call__(self, plot):
        plot._axes.hold(True)
        ax = plot.data.axis
        xax = plot.data.ds.coordinates.x_axis[ax]
        yax = plot.data.ds.coordinates.y_axis[ax]

        if not hasattr(self.vertices, "in_units"):
            vertices = plot.data.pf.arr(self.vertices, "code_length")
        else:
            vertices = self.vertices
        l_cy = triangle_plane_intersect(plot.data.axis, plot.data.coord, vertices)[:,:,(xax, yax)]
        # reformat for conversion to plot coordinates
        l_cy = np.rollaxis(l_cy,0,3)
        # convert all line starting points
        l_cy[0] = self.convert_to_plot(plot,l_cy[0])
        l_cy[1] = self.convert_to_plot(plot,l_cy[1])
        # convert all line ending points
        l_cy = np.rollaxis(l_cy,2,0)
        # create line collection and add it to the plot
        lc = matplotlib.collections.LineCollection(l_cy, **self.plot_args)
        plot._axes.add_collection(lc)
        plot._axes.hold(False)

class TimestampCallback(PlotCallback):
    """
    annotate_timestamp(corner='lowerleft', time=True, redshift=False, 
                       time_format="t = {time:.0f} {units}", time_unit=None, 
                       redshift_format="z = {redshift:.2f}", 
                       bbox=False, pos=None, text_args=None, bbox_args=None)

    Annotates the timestamp and/or redshift of the data output at a specified
    location in the image (either in a present corner, or by specifying (x,y)
    image coordinates with the pos argument.  If no time_units are specified, 
    it will automatically choose appropriate units.  It allows for custom 
    formatting of the time and redshift information, as well as the 
    specification of a bounding box around the text.

    Parameters
    ----------
    corner : string, optional
        Corner sets up one of 4 predeterimined locations for the timestamp
        to be displayed in the image: 'upperleft', 'upperright', 'lowerleft',
        'lowerright' (also allows None). This value will be trumped by the 
        optional 'pos' keyword.
    time : boolean, optional
        Whether or not to show the ds.current_time of the data output.  Can
        be used solo or in conjunction with redshift parameter.
    redshift : boolean, optional
        Whether or not to show the ds.current_time of the data output.  Can
        be used solo or in conjunction with time parameter.
    time_format : string, optional
        This specifies the format of the time output assuming "time" is the 
        number of time and "unit" is units of the time (e.g. 's', 'Myr', etc.)
        The time can be specified to arbitrary precision according to printf
        formatting codes (defaults to .1f -- a float with 1 digits after 
        decimal).
    time_unit : string, optional
        time_unit must be a valid yt time unit (e.g. 's', 'min', 'hr', 'yr', 
        'Myr', etc.)
    redshift_format : string, optional
        This specifies the format of the redshift output.  The redshift can
        be specified to arbitrary precision according to printf formatting 
        codes (defaults to 0.2f -- a float with 2 digits after decimal).
    bbox : boolean, optional
        Whether or not a bounding box should be included around the text
        If so, it uses the bbox_args to set the matplotlib FancyBboxPatch 
        object.  
    pos : tuple of floats, optional
        The image location of the timestamp in image coords (i.e. (x,y) = 
        (0..1, 0..1).  Setting pos trumps the corner parameter.
    text_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        text object.  Defaults: {'color':'white'}.
    bbox_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        FancyBboxPatch object as the bounding box around the text.  
        Defaults: {'boxstyle':'square,pad=0.3', 'facecolor':'black', 
                  'linewidth':3, 'edgecolor':'white', 'alpha':'0.5'}
    """
    _type_name = "timestamp"
    # Defaults
    _text_args = {'color':'white'}
    _bbox_args = {'boxstyle':'square,pad=0.3', 'facecolor':'black', 
                  'linewidth':3, 'edgecolor':'white', 'alpha':0.5}

    def __init__(self, corner='lowerleft', time=True, redshift=False, 
                 time_format="t = {time:.1f} {units}", time_unit=None,
                 redshift_format="z = {redshift:.2f}", bbox=False,
                 pos=None, text_args=None, bbox_args=None):

        # Set position based on corner argument.
        self.corner = corner
        self.time = time
        self.redshift = redshift
        self.time_format = time_format
        self.redshift_format = redshift_format
        self.time_unit = time_unit
        self.pos = pos
        if text_args is None: self.text_args = self._text_args
        else: self.text_args = text_args
        if bbox_args is None: self.bbox_args = self._bbox_args
        else: self.bbox_args = bbox_args

        # if bbox is not desired, set bbox_args to {}
        if not bbox: self.bbox_args = {}

    def __call__(self, plot):
        # Setting pos trumps corner argument
        if self.pos is None:
            if self.corner == 'upperleft':
                self.pos = (0.03, 0.97)
                self.text_args['horizontalalignment'] = 'left'
                self.text_args['verticalalignment'] = 'top'
            elif self.corner == 'upperright':
                self.pos = (0.97, 0.97)
                self.text_args['horizontalalignment'] = 'right'
                self.text_args['verticalalignment'] = 'top'
            elif self.corner == 'lowerleft':
                self.pos = (0.03, 0.03)
                self.text_args['horizontalalignment'] = 'left'
                self.text_args['verticalalignment'] = 'bottom'
            elif self.corner == 'lowerright':
                self.pos = (0.97, 0.03)
                self.text_args['horizontalalignment'] = 'right'
                self.text_args['verticalalignment'] = 'bottom'
            elif self.corner is None:
                self.pos = (0.5, 0.5)
                self.text_args['horizontalalignment'] = 'center'
                self.text_args['verticalalignment'] = 'center'
            else:
                print "Argument 'corner' must be set to 'upperleft'," \
                      "'upperright', 'lowerleft', 'lowerright', or None"
                raise SyntaxError

        self.text = ""

        # If we're annotating the time, put it in the correct format
        if self.time:

            # If no time_units are set, then identify a best fit time unit
            if self.time_unit is None:
                t = plot.data.ds.current_time.in_units('s')
                scale_keys = ['fs', 'ps', 'ns', 'us', 'ms', 's', 'hr', 
                              'day', 'yr', 'kyr', 'Myr', 'Gyr', 'Tyr', 
                              'Pyr', 'Eyr', 'Zyr', 'Yyr']
                for i, k in enumerate(scale_keys):
                    if t < YTQuantity(1, k):
                        break
                    t.convert_to_units(k)
                self.time_unit = scale_keys[i-1]
            # If time_unit is set, use it
            else:
                t = plot.data.ds.current_time.in_units(self.time_unit)
            self.text += self.time_format.format(time=float(t), 
                                                 units=self.time_unit)

        # If time and redshift both shown, do one on top of the other
        if self.time and self.redshift:
            self.text += "\n"

        # If we're annotating the redshift, put it in the correct format
        if self.redshift:
            try:
                z = np.abs(plot.data.ds.current_redshift)
            except AttributeError:
                print "Dataset does not have current_redshift. Set redshift=False."
                raise AttributeError
            self.text += self.redshift_format.format(redshift=float(z))

        # This is just a fancy wrapper around the TextLabelCallback
        tcb = TextLabelCallback(self.pos, self.text, 
                                text_args=self.text_args, 
                                bbox_args=self.bbox_args)
        return tcb(plot)

class ScaleCallback(PlotCallback):
    """
    annotate_scale(corner='lowerright', coeff=None, unit=None, pos=None,
                   max_frac=0.2, min_frac=0.018,
                   text_args=None, plot_args=None)

    Annotates the scale of the plot at a specified location in the image
    (either in a preset corner, or by specifying (x,y) image coordinates with
    the pos argument.  Coeff and units (e.g. 1 Mpc) refer to the distance scale
    you desire to show on the plot.  If no coeff and units are specified, 
    an appropriate pair will be determined such that your scale bar is never
    smaller than min_frac or greater than max_frac of your plottable axis 
    length.  For additional text and plot arguments for the text and line,
    include them as dictionaries to pass to text_args and plot_args.
    
    Parameters
    ----------
    corner : string, optional
        Corner sets up one of 4 predeterimined locations for the timestamp
        to be displayed in the image: 'upperleft', 'upperright', 'lowerleft',
        'lowerright' (also allows None). This value will be trumped by the 
        optional 'pos' keyword.
    coeff : float, optional
        The coefficient before the unit scale.  If set to None along with
        unit keyword, will be automatically determined to be a power of 10
        relative to the best-fit scale.
    unit : string, optional
        unit must be a valid yt distance unit (e.g. 'm', 'km', 'AU', 'pc', 
        'kpc', etc.) or set to None.  If set to None, will be automatically
        determined to be the best-fit to the data.
    pos : tuple of floats, optional
        The image location of the timestamp in image coords (i.e. (x,y) = 
        (0..1, 0..1).  Setting pos trumps the corner parameter.
    min_frac, max_frac: float, optional
        The minimum/maximum fraction of the axis width for the scale bar to 
        extend. A value of 1 would allow the scale bar to extend across the
        entire axis width.  Only used for automatically calculating 
        best-fit coeff and unit when neither is specified, otherwise 
        disregarded.
    text_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        text object.  Defaults: {'color':'white', 
        'horizontalalignment':'center', 'verticalalignment':'top'}.
    plot_args : dictionary, optional
        A dictionary of any arbitrary parameters to be passed to the Matplotlib
        line object.  Defaults: {'color':'white', 'linewidth':3}.
    """
    _type_name = "scale"
    # Defaults
    _text_args = {'horizontalalignment':'center', \
                  'verticalalignment':'top', \
                  'color':'white'}
    _plot_args = {'color':'white', 'linewidth':3}

    def __init__(self, corner='lowerright', coeff=None, unit=None, pos=None, 
                 max_frac=0.20, min_frac=0.018,
                 text_args=None, plot_args=None):

        # Set position based on corner argument.
        self.corner = corner
        self.coeff = coeff
        self.unit = unit
        self.pos = pos
        self.max_frac = max_frac
        self.min_frac = min_frac
        if text_args is None: self.text_args = self._text_args
        else: self.text_args = text_args
        if plot_args is None: self.plot_args = self._plot_args
        else: self.plot_args = plot_args

    def __call__(self, plot):
        # Callback only works for plots with axis ratios of 1
        xsize = plot.xlim[1] - plot.xlim[0]
        ysize = plot.ylim[1] - plot.ylim[0]
        if xsize != ysize:
            raise RuntimeError("Scale callback only works for plots with", \
                                "axis ratios of 1: xsize = %s, ysize = %s." %
                                (xsize, ysize))

        # Setting pos trumps corner argument
        if self.pos is None:
            if self.corner == 'upperleft':
                self.pos = (0.12, 0.971)
            elif self.corner == 'upperright':
                self.pos = (0.88, 0.971)
            elif self.corner == 'lowerleft':
                self.pos = (0.12, 0.062)
            elif self.corner == 'lowerright':
                self.pos = (0.88, 0.062)
            elif self.corner is None:
                self.pos = (0.5, 0.5)
            else:
                print "Argument 'corner' must be set to 'upperleft',", \
                      "'upperright', 'lowerleft', 'lowerright', or None"
                raise SyntaxError

        # When identifying a best fit distance unit, do not allow scale marker
        # to be greater than max_frac fraction of xaxis or under min_frac 
        # fraction of xaxis 
        max_scale = self.max_frac * xsize
        min_scale = self.min_frac * xsize

        if self.coeff is None:
            self.coeff = 1.

        # If no units are set, then identify a best fit distance unit
        if self.unit is None:
            scale_keys = ['fm', 'pm', 'nm', 'um', 'mm', 'cm', 'm', 
                          'km', 'AU', 'pc', 'kpc', 'Mpc', 
                          'Gpc', 'Tpc', 'Ppc', 'Epc', 'Zpc', 'Ypc']
            for i, k in enumerate(scale_keys):
                if i == len(scale_keys)-1: break
                # How many OOM between this scale and next? Use each as 
                # coefficient
                powers = np.ceil(np.log10(YTQuantity(1, scale_keys[i+1]) / \
                                 YTQuantity(1, scale_keys[i])))
                coefficients = 10**np.arange(powers)
                for j in coefficients:
                    if not (min_scale < YTQuantity(j, k) and max_scale > YTQuantity(j, k)):
                        continue
                    self.coeff = j
                    self.unit = scale_keys[i]
                    break
        self.scale = YTQuantity(self.coeff, self.unit)
        self.text = "{scale} {units}".format(scale=int(self.coeff), 
                                             units=self.unit)
        image_scale = (plot.frb.convert_distance_x(self.scale) / \
                       plot.frb.convert_distance_x(xsize)).v

        # This is just a fancy wrapper around the TextLabelCallback and the
        # ImageLineCallback
        pos_line_start = (self.pos[0]-image_scale/2, self.pos[1]+0.01)
        pos_line_end = (self.pos[0]+image_scale/2, self.pos[1]+0.01)
        icb = ImageLineCallback(pos_line_start, pos_line_end, plot_args=self.plot_args)
        icb(plot)
        tcb = TextLabelCallback(self.pos, self.text, 
                                text_args=self.text_args)
        return tcb(plot)

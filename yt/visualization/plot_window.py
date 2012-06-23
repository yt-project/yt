"""
A plotting mechanism based on the idea of a "window" into the data.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 J. S. Oishi.  All Rights Reserved.

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
import base64
import matplotlib.pyplot
import cStringIO
from functools import wraps

import numpy as na
from .image_writer import \
    write_image, apply_colormap
from .fixed_resolution import \
    FixedResolutionBuffer
from .plot_modifications import get_smallest_appropriate_unit
from .tick_locators import LogLocator, LinearLocator

from yt.funcs import *
from yt.utilities.lib import write_png_to_string
from yt.utilities.definitions import \
    x_dict, x_names, \
    y_dict, y_names, \
    axis_names, \
    axis_labels

def invalidate_data(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._data_valid = False
        args[0]._plot_valid = False
        args[0]._recreate_frb()
        if args[0]._initfinished:
            args[0]._setup_plots()
        return rv
    return newfunc

def invalidate_plot(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._plot_valid = False
        args[0]._setup_plots()
        return rv
    return newfunc

field_transforms = {}

class FieldTransform(object):
    def __init__(self, name, func, locator):
        self.name = name
        self.func = func
        self.locator = locator
        field_transforms[name] = self

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def ticks(self, mi, ma):
        try:
            ticks = self.locator(mi, ma)
        except:
            ticks = []
        return ticks

log_transform = FieldTransform('log10', na.log10, LogLocator())
linear_transform = FieldTransform('linear', lambda x: x, LinearLocator())

def SlicePlot(pf, axis, fields, center=None, width=None, origin='center-window'):
    (bounds,center) = GetBoundsAndCenter(axis,center,width,pf)
    slice = pf.h.slice(axis,center[axis],fields=fields)
    return PWViewerMPL(slice,bounds,origin=origin)

def ProjectionPlot(pf, axis, fields, center=None, width=None,
                   weight_field=None, max_level=None, origin='center-window'):
    (bounds,center) = GetBoundsAndCenter(axis,center,width,pf)
    proj = pf.h.proj(axis,fields,weight_field=weight_field,max_level=max_level,center=center)
    return PWViewerMPL(proj,bounds,origin=origin)

def GetBoundsAndCenter(axis, center, width, pf):
    if width == None:
        width = (pf.domain_right_edge - pf.domain_left_edge)
    if iterable(width):
        w,u = width
        width = w/pf[u]
    if center == None:
        v, center = pf.h.mind_max("Density")
    elif center == "center" or center == "c":
        center = (pf.domain_right_edge + pf.domain_left_edge)/2.0
    bounds = [center[x_dict[axis]]-width/2,
              center[x_dict[axis]]+width/2,
              center[y_dict[axis]]-width/2,
              center[y_dict[axis]]+width/2] 
    return (bounds,center)

class PlotWindow(object):
    _plot_valid = False
    _colorbar_valid = False
    _contour_info = None
    _vector_info = None
    def __init__(self, data_source, bounds, buff_size=(800,800), antialias = True, 
                 periodic = True, origin='center-window'):
        r"""
        PlotWindow(data_source, bounds, buff_size=(800,800), antialias = True)
        
        A ploting mechanism based around the concept of a window into a
        data source. It can have arbitrary fields, each of which will be
        centered on the same viewpoint, but will have individual zlimits. 
        
        The data and plot are updated separately, and each can be
        invalidated as the object is modified.
        
        Data is handled by a FixedResolutionBuffer object.

        Parameters
        ----------
        data_source : :class:`yt.data_objects.data_containers.AMRProjBase` or :class:`yt.data_objects.data_containers.AMRSliceBase`
            This is the source to be pixelized, which can be a projection or a
            slice.  (For cutting planes, see
            `yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer`.)
        bounds : sequence of floats
            Bounds are the min and max in the image plane that we want our
            image to cover.  It's in the order of (xmin, xmax, ymin, ymax),
            where the coordinates are all in the appropriate code units.
        buff_size : sequence of ints
            The size of the image to generate.
        antialias : boolean
            This can be true or false.  It determines whether or not sub-pixel
            rendering is used during data deposition.

        """
        self._initfinished = False
        self.center = None
        self.plots = {}
        self._periodic = periodic
        self.data_source = data_source
        self.buff_size = buff_size
        self.antialias = True
        self.set_window(bounds) # this automatically updates the data and plot
        self.origin = origin
        if self.data_source.center is not None:
            center = [self.data_source.center[i] for i in range(len(self.data_source.center)) if i != self.data_source.axis]
            self.set_center(center)
        self._initfinished = True

    def __getitem__(self, item):
        return self.plots[item]

    def _recreate_frb(self):
        try:
            bounds = self.bounds
            self._frb = FixedResolutionBuffer(self.data_source, 
                                              bounds, self.buff_size, 
                                              self.antialias, periodic=self._periodic)
        except:
            raise RuntimeError("Failed to repixelize.")
        self._frb._get_data_source_fields()
        self.pf = self._frb.pf
        self._data_valid = True
        
    def _setup_plots(self):
        pass

    @property
    def fields(self):
        return self._frb.data.keys()

    @property
    def width(self):
        Wx = self.xlim[1] - self.xlim[0]
        Wy = self.ylim[1] - self.ylim[0]
        return (Wx, Wy)

    @property
    def bounds(self):
        return self.xlim+self.ylim

    @invalidate_data
    def zoom(self, factor):
        r"""This zooms the window by *factor*.

        Parameters
        ----------
        factor : float
            multiplier for the current width

        """
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        nWx, nWy = Wx/factor, Wy/factor
        self.xlim = (centerx - nWx*0.5, centerx + nWx*0.5)
        self.ylim = (centery - nWy*0.5, centery + nWy*0.5)

    @invalidate_data
    def pan(self, deltas):
        r"""Pan the image by specifying absolute code unit coordinate deltas.
        
        Parameters
        ----------
        deltas : sequence of floats
            (delta_x, delta_y) in *absolute* code unit coordinates

        """
        self.xlim = (self.xlim[0] + deltas[0], self.xlim[1] + deltas[0])
        self.ylim = (self.ylim[0] + deltas[1], self.ylim[1] + deltas[1])

    @invalidate_data
    def pan_rel(self, deltas):
        r"""Pan the image by specifying relative deltas, to the FOV.
        
        Parameters
        ----------
        deltas : sequence of floats
            (delta_x, delta_y) in *relative* code unit coordinates

        """
        Wx, Wy = self.width
        self.xlim = (self.xlim[0] + Wx*deltas[0], self.xlim[1] + Wx*deltas[0])
        self.ylim = (self.ylim[0] + Wy*deltas[1], self.ylim[1] + Wy*deltas[1])

    @invalidate_data
    def set_field(self):
        pass

    @invalidate_data
    def set_window(self, bounds):
        if self.center is not None:
            dx = bounds[1] - bounds[0]
            dy = bounds[3] - bounds[2]
            self.xlim = (self.center[0] - dx/2., self.center[0] + dx/2.)
            self.ylim = (self.center[1] - dy/2., self.center[1] + dy/2.)
            mylog.info("xlim = %f %f" %self.xlim)
            mylog.info("ylim = %f %f" %self.ylim)
        else:
            self.xlim = bounds[0:2]
            self.ylim = bounds[2:]
            
    @invalidate_data
    def set_width(self, new_width):
        """set the width of the plot window

        parameters
        ----------
        new_width : float
            the width of the image in code units.

        """
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        self.xlim = (centerx - new_width/2.,
                     centerx + new_width/2.)
        self.ylim = (centery - new_width/2.,
                     centery + new_width/2.)

    @invalidate_data
    def set_center(self, new_center):
        if new_center is None:
            self.center = None
        else:
            self.center = new_center
        self.set_window(self.bounds)

    @property
    def width(self):
        Wx = self.xlim[1] - self.xlim[0]
        Wy = self.ylim[1] - self.ylim[0]
        return (Wx, Wy)

    # @property
    # def window(self):
    #     return self.xlim + self.ylim

    @invalidate_data
    def set_antialias(self,aa):
        self.antialias = aa

    @invalidate_plot
    def set_contour_info(self, field_name, n_cont = 8, colors = None,
                         logit = True):
        if field_name == "None" or n_cont == 0:
            self._contour_info = None
            return
        self._contour_info = (field_name, n_cont, colors, logit)

    @invalidate_plot
    def set_vector_info(self, skip, scale = 1):
        self._vector_info = (skip, scale)

class PWViewer(PlotWindow):
    """A viewer for PlotWindows.

    """
    def __init__(self, *args,**kwargs):
        setup = kwargs.pop("setup", True)
        PlotWindow.__init__(self, *args,**kwargs)
        self._field_transform = {}
        self._colormaps = defaultdict(lambda: 'algae')
        for field in self._frb.data.keys():
            if self.pf.field_info[field].take_log:
                self._field_transform[field] = log_transform
            else:
                self._field_transform[field] = linear_transform

        if setup: self._setup_plots()

    @invalidate_plot
    def set_log(self, field, log):
        """set a field to log or linear.
        
        Parameters
        ----------
        field : string
            the field to set a transform
        log : boolean
            Log on/off.

        """
        if log:
            self._field_transform[field] = log_transform
        else:
            self._field_transform[field] = linear_transform

    @invalidate_plot
    def set_transform(self, field, name):
        if name not in field_transforms: 
            raise KeyError(name)
        self._field_transform[field] = field_transforms[name]

    @invalidate_plot
    def set_cmap(self, field, cmap_name):
        self._colorbar_valid = False
        self._colormaps[field] = cmap_name

    @invalidate_plot
    def set_zlim(self):
        pass

    def get_metadata(self, field, strip_mathml = True, return_string = True):
        fval = self._frb[field]
        mi = fval.min()
        ma = fval.max()
        x_width = self.xlim[1] - self.xlim[0]
        y_width = self.ylim[1] - self.ylim[0]
        unit = get_smallest_appropriate_unit(x_width, self.pf)
        units = self.get_field_units(field, strip_mathml)
        center = getattr(self._frb.data_source, "center", None)
        if center is None or self._frb.axis == 4:
            xc, yc, zc = -999, -999, -999
        else:
            center[x_dict[self._frb.axis]] = 0.5 * (
                self.xlim[0] + self.xlim[1])
            center[y_dict[self._frb.axis]] = 0.5 * (
                self.ylim[0] + self.ylim[1])
            xc, yc, zc = center
        if return_string:
            md = _metadata_template % dict(
                pf = self.pf,
                x_width = x_width*self.pf[unit],
                y_width = y_width*self.pf[unit],
                unit = unit, units = units, mi = mi, ma = ma,
                xc = xc, yc = yc, zc = zc)
        else:
            md = dict(pf = self.pf,
                      x_width = x_width*self.pf[unit],
                      y_width = y_width*self.pf[unit],
                      unit = unit, units = units, mi = mi, ma = ma,
                      xc = xc, yc = yc, zc = zc)
        return md

    def get_field_units(self, field, strip_mathml = True):
        ds = self._frb.data_source
        pf = self.pf
        if ds._type_name == "slice":
            units = pf.field_info[field].get_units()
        elif ds._type_name == "proj":
            units = pf.field_info[field].get_projected_units()
        else:
            units = ""
        if strip_mathml:
            units = units.replace(r"\rm{", "").replace("}","")
        return units


class PWViewerMPL(PWViewer):
    """Viewer using matplotlib as a backend via the WindowPlotMPL. 

    """
    _current_field = None

    def _setup_plots(self):
        if self._current_field is not None:
            fields = [self._current_field]
        else:
            fields = self._frb.data.keys()
        self._colorbar_valid = True
        for f in self.fields:
            md = self.get_metadata(f, strip_mathml = False, return_string = False)
            if self.origin == 'center-window':
                extent = [self.xlim[i] - (self.xlim[0]+self.xlim[1])/2. for i in (0,1)]
                extent.extend([self.ylim[i] - (self.ylim[0]+self.ylim[1])/2. for i in (0,1)])
            elif self.origin == 'center-domain':
                pass
            elif self.origin == 'left-domain':
                pass
            else:
                raise RuntimeError('Origin keyword not recognized')

            extent = [el*self.pf[md['unit']] for el in extent]

            self.plots[f] = WindowPlotMPL(self._frb[f], extent)
            
            cb = matplotlib.pyplot.colorbar(self.plots[f].image,cax = self.plots[f].cax)
            axis_index = self.data_source.axis
            
            labels = [r'$\rm{'+axis_labels[axis_index][i].encode('string-escape')+
                      r'\/\/('+md['unit'].encode('string-escape')+r')}$' for i in (0,1)]

            self.plots[f].axes.set_xlabel(labels[0])
            self.plots[f].axes.set_ylabel(labels[1])
             
            cb.set_label(r'$\rm{'+f.encode('string-escape')+r'}\/\/('+md['units']+r')$')

        self._plot_valid = True

    def save(self,name):
        for k,v in self.plots.iteritems():
            n = "%s_%s" % (name, k)
            v.save(n)

class PWViewerRaw(PWViewer):
    """A PlotWindow viewer that writes raw pngs (no MPL, no axes).

    """
    def _setup_plots(self):
        self.save('')
        self._plot_valid = True

    def save(self,name):
        for field in self._frb.data.keys():
            nm = "%s_%s.png" % (name,field)
            print "writing %s" % nm
            write_image(self._frb[field],nm)

_metadata_template = """
%(pf)s<br>
<br>
Field of View:  %(x_width)0.3f %(unit)s<br>
Minimum Value:  %(mi)0.3e %(units)s<br>
Maximum Value:  %(ma)0.3e %(units)s<br>
Central Point:  (data coords)<br>
&nbsp;&nbsp;&nbsp;%(xc)0.14f<br>
&nbsp;&nbsp;&nbsp;%(yc)0.14f<br>
&nbsp;&nbsp;&nbsp;%(zc)0.14f
"""

class PWViewerExtJS(PWViewer):
    """A viewer for the web interface.

    """
    _ext_widget_id = None
    _current_field = None
    _widget_name = "plot_window"

    def _setup_plots(self):
        from yt.gui.reason.bottle_mods import PayloadHandler
        ph = PayloadHandler()
        if self._current_field is not None \
           and self._ext_widget_id is not None:
            fields = [self._current_field]
            addl_keys = {'type': 'widget_payload',
                         'widget_id': self._ext_widget_id}
        else:
            fields = self._frb.data.keys()
            addl_keys = {}
        if self._colorbar_valid == False:
            addl_keys['colorbar_image'] = self._get_cbar_image()
            self._colorbar_valid = True
        min_zoom = 200*self.pf.h.get_smallest_dx() * self.pf['unitary']
        for field in fields:
            to_plot = apply_colormap(self._frb[field],
                func = self._field_transform[field],
                cmap_name = self._colormaps[field])
            pngs = self._apply_modifications(to_plot)
            img_data = base64.b64encode(pngs)
            # We scale the width between 200*min_dx and 1.0
            x_width = self.xlim[1] - self.xlim[0]
            zoom_fac = na.log10(x_width*self.pf['unitary'])/na.log10(min_zoom)
            zoom_fac = 100.0*max(0.0, zoom_fac)
            ticks = self.get_ticks(field)
            payload = {'type':'png_string',
                       'image_data':img_data,
                       'metadata_string': self.get_metadata(field),
                       'zoom': zoom_fac,
                       'ticks': ticks}
            payload.update(addl_keys)
            ph.add_payload(payload)

    def _apply_modifications(self, img):
        if self._contour_info is None and self._vector_info is None:
            return write_png_to_string(img)
        from matplotlib.figure import Figure
        from yt.visualization._mpl_imports import \
            FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS
        from yt.utilities.delaunay.triangulate import Triangulation as triang

        vi, vj, vn = img.shape

        # Now we need to get our field values
        fig = Figure((vi/100.0, vj/100.0), dpi = 100)
        fig.figimage(img)
        # Add our contour
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
        ax.patch.set_alpha(0.0)

        # Now apply our modifications
        self._apply_contours(ax, vi, vj)
        self._apply_vectors(ax, vi, vj)

        canvas = FigureCanvasAgg(fig)
        f = cStringIO.StringIO()
        canvas.print_figure(f)
        f.seek(0)
        img = f.read()
        return img

    def _apply_contours(self, ax, vi, vj):
        if self._contour_info is None: return 
        plot_args = {}
        field, number, colors, logit = self._contour_info
        if colors is not None: plot_args['colors'] = colors

        raw_data = self._frb.data_source
        b = self._frb.bounds
        xi, yi = na.mgrid[b[0]:b[1]:(vi / 8) * 1j,
                          b[2]:b[3]:(vj / 8) * 1j]
        x = raw_data['px']
        y = raw_data['py']
        z = raw_data[field]
        if logit: z = na.log10(z)
        fvals = triang(x,y).nn_interpolator(z)(xi,yi).transpose()[::-1,:]

        ax.contour(fvals, number, colors='w')
        
    def _apply_vectors(self, ax, vi, vj):
        if self._vector_info is None: return 
        skip, scale = self._vector_info

        nx = self._frb.buff_size[0]/skip
        ny = self._frb.buff_size[1]/skip
        new_frb = FixedResolutionBuffer(self._frb.data_source,
                        self._frb.bounds, (nx,ny))

        axis = self._frb.data_source.axis
        fx = "%s-velocity" % (axis_names[x_dict[axis]])
        fy = "%s-velocity" % (axis_names[y_dict[axis]])
        px = new_frb[fx][::-1,:]
        py = new_frb[fy][::-1,:]
        x = na.mgrid[0:vi-1:ny*1j]
        y = na.mgrid[0:vj-1:nx*1j]
        # Always normalize, then we scale
        nn = ((px**2.0 + py**2.0)**0.5).max()
        px /= nn
        py /= nn
        print scale, px.min(), px.max(), py.min(), py.max()
        ax.quiver(x, y, px, py, scale=float(vi)/skip)
        
    def get_ticks(self, field, height = 400):
        # This will eventually change to work with non-logged fields
        ticks = []
        transform = self._field_transform[field]
        mi, ma = self._frb[field].min(), self._frb[field].max()
        tick_locs = transform.ticks(mi, ma)
        mi, ma = transform((mi, ma))
        for v1,v2 in zip(tick_locs, transform(tick_locs)):
            if v2 < mi or v2 > ma: continue
            p = height - height * (v2 - mi)/(ma - mi)
            ticks.append((p,v1,v2))
        return ticks

    def _get_cbar_image(self, height = 400, width = 40, field = None):
        if field is None: field = self._current_field
        cmap_name = self._colormaps[field]
        vals = na.mgrid[1:0:height * 1j] * na.ones(width)[:,None]
        vals = vals.transpose()
        to_plot = apply_colormap(vals, cmap_name = cmap_name)
        pngs = write_png_to_string(to_plot)
        img_data = base64.b64encode(pngs)
        return img_data

    # This calls an invalidation routine from within
    def scroll_zoom(self, value):
        # We accept value from 0..100, and assume it has been set from the
        # scroll bar.  In that case, we undo the logic for calcualting
        # 'zoom_fac' from above.
        min_val = 200*self.pf.h.get_smallest_dx()
        unit = self.pf['unitary']
        width = (min_val**(value/100.0))/unit
        self.set_width(width)

    def image_recenter(self, img_x, img_y, img_size_x, img_size_y):
        dx = (self.xlim[1] - self.xlim[0]) / img_size_x
        dy = (self.ylim[1] - self.ylim[0]) / img_size_y
        new_x = img_x * dx + self.xlim[0]
        new_y = img_y * dy + self.ylim[0]
        print img_x, img_y, dx, dy, new_x, new_y
        self.set_center((new_x, new_y))

    @invalidate_data
    def set_current_field(self, field):
        self._current_field = field
        self._frb[field]
        if self.pf.field_info[field].take_log:
            self._field_transform[field] = log_transform
        else:
            self._field_transform[field] = linear_transform

class PlotMPL(object):
    """A base class for all yt plots made using matplotlib.

    YtPlot and the classes that derive from it are *by design* limited
    and designed for rapid, production quality plot production, rather
    than full generality. If you require more customization of the end
    result, these objects are designed to return to you the basic data
    so you the user can insert them into a matplotlib figure on your
    own outside of the YtPlot class.

    """
    datalabel = None
    figure = None
    def __init__(self, field, size):
        self._plot_valid = True
        self.figure = matplotlib.pyplot.figure(figsize=size,frameon=True)
        self.axes = self.figure.add_axes((.07,.10,.8,.8))
        self.cax = self.figure.add_axes((.86,.10,.04,.8))

    def save(self,name):
        print "saving plot %s" % name
        self.figure.savefig('%s.png' % name)

class WindowPlotMPL(PlotMPL):
    zmin = None
    zmax = None
    zlabel = None

    def __init__(self, data, extent, size=(9,8)):
        PlotMPL.__init__(self, data, size)
        self.__init_image(data, extent)

    def __init_image(self, data, extent):
        self.image = self.axes.imshow(data,origin='lower',extent=extent,
                                      norm=matplotlib.colors.LogNorm())

    @invalidate_plot
    def set_zlim(self, zmin, zmax):
        self.zmin = zmin
        self.zmax = zmax


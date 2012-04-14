"""
This is a simple mechanism for interfacing with Profile and Phase plots

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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
import types

from functools import wraps
import numpy as na

from .image_writer import \
    write_image, apply_colormap
from yt.utilities.amr_utils import \
    write_png_to_string
from yt.data_objects.profiles import \
    BinnedProfile1D, \
    BinnedProfile2D
from .plot_types import ProfilePlot, PhasePlot
from .tick_locators import LogLocator, LinearLocator
from yt.utilities.logger import ytLogger as mylog

def invalidate_plot(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        args[0]._plot_valid = False
        args[0]._setup_plot()
        return f(*args, **kwargs)
    return newfunc

class AxisSpec(object):
    title = None
    bounds = None
    scale = None
    ticks = None

    def calculate_ticks(self):
        if self.scale == 'log':
            locator = LogLocator()
        elif self.scale == 'linear':
            locator = LinearLocator()
        else:
            raise NotImplementedError
        self.ticks = locator(*self.bounds)

class ColorbarSpec(AxisSpec):
    cmap = None
    display = True

class PlotContainer(object):
    x_spec = None
    y_spec = None
    x_values = None
    y_values = None
    plot_spec = None

    def to_mpl(self, place = None):
        import _mpl_imports as mpl
        if isinstance(place, mpl.matplotlib.figure.Figure):
            figure, place = place, None
            place = None
        else:
            figure = mpl.matplotlib.figure.Figure((10,8))
        if isinstance(place, mpl.matplotlib.axes.Axes):
            axes, place = place, None
        else:
            axes = figure.add_subplot(1,1,1)
        if self.x_spec.scale == 'log' and \
           self.y_spec.scale == 'log':
            func = axes.loglog
        elif self.x_spec.scale == 'log':
            func = axes.semilogx
        elif self.y_spec.scale == 'log':
            func = axes.semilogy
        else:
            func = axes.plot
        if self.plot_spec is None:
            kwargs = {}
        else:
            kwargs = self.plot_spec
        func(self.x_values, self.y_values, **kwargs)
        if self.x_spec.title is not None:
            axes.set_xlabel(self.x_spec.title)
        if self.y_spec.title is not None:
            axes.set_ylabel(self.y_spec.title)
        if isinstance(place, types.StringTypes):
            canvas = mpl.FigureCanvasAgg(figure)
            canvas.print_figure(place)
        return figure, axes

class ImagePlotContainer(object):
    x_spec = None
    y_spec = None
    image = None
    cbar = None

    def to_mpl(self, place = None):
        import _mpl_imports as mpl
        if isinstance(place, mpl.matplotlib.figure.Figure):
            figure, place = place, None
            place = None
        else:
            figure = mpl.matplotlib.figure.Figure((10,8))
        if isinstance(place, mpl.matplotlib.axes.Axes):
            axes, place = place, None
        else:
            axes = figure.add_subplot(1,1,1)
        # We'll go with a mesh here, even if it's inappropriate
        use_mesh = False
        xmi, xma = self.x_spec.bounds
        if self.x_spec.scale == 'log':
            x_bins = na.logspace(na.log10(xmi), na.log10(xma),
                                 self.image.shape[0]+1)
            use_mesh = True
        else:
            x_bins = na.logspace(xmi, xma, self.image.shape[0]+1)

        ymi, yma = self.y_spec.bounds
        if self.y_spec.scale == 'log':
            y_bins = na.logspace(na.log10(ymi), na.log10(yma),
                                 self.image.shape[0]+1)
            use_mesh = True
        else:
            y_bins = na.logspace(ymi, yma, self.image.shape[0]+1)

        im = self.image
        if self.cbar.scale == 'log':
            norm = mpl.matplotlib.colors.LogNorm()
        else:
            norm = mpl.matplotlib.colors.Normalize()
        if use_mesh:
            pcm = axes.pcolormesh(x_bins, y_bins, self.image, norm=norm,
                                  shading='flat', cmap = self.cbar.cmap,
                                  rasterized=True)
            if self.x_spec.scale == 'log': axes.set_xscale("log")
            if self.y_spec.scale == 'log': axes.set_yscale("log")
        else:
            axes.imshow(self.image, origin='lower', interpolation='nearest',
                        cmap = self.cbar.cmap, extent = [xmi,xma,ymi,yma],
                        norm = norm)
        if self.x_spec.title is not None:
            axes.set_xlabel(self.x_spec.title)
        if self.y_spec.title is not None:
            axes.set_ylabel(self.y_spec.title)
        if isinstance(place, types.StringTypes):
            canvas = mpl.FigureCanvasAgg(figure)
            canvas.print_figure(place)
        return figure, axes

class PhasePlotter(object):
    scale = None
    _current_field = None

    def __init__(self, data_source, field_x, field_y, field_z,
                 weight="CellMassMsun", accumulation=False,
                 x_bins=128, x_log=True, x_bounds=None,
                 y_bins=128, y_log=True, y_bounds=None,
                 lazy_reader=True, fractional=False):
        r"""From an existing object, create a 2D, binned profile.

        This function will accept an existing `AMRData` source and from that,
        it will generate a `Binned2DProfile`, based on the specified options.
        This is useful if you have extracted a region, or if you wish to bin
        some set of massages data -- or even if you wish to bin anything other
        than a sphere.  The profile will be 2D, which means while it can have
        an arbitrary number of fields, those fields will all be binned based on
        two fields.

        Parameters
        ----------
        data_source : `yt.data_objects.api.AMRData`
            This is a data source respecting the `AMRData` protocol (i.e., it
            has grids and so forth) that will be used as input to the profile
            generation.
        fields : list of strings
            The first element of this list is the field by which we will bin
            into the x-axis, the second is the field by which we will bin onto
            the y-axis.  All subsequent fields will be binned and their
            profiles added to the underlying `BinnedProfile2D`.
        cmap : string, optional
            An acceptable colormap.  See either yt.visualization.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        weight : string, default "CellMassMsun"
            The weighting field for an average.  This defaults to mass-weighted
            averaging.
        accumulation : list of booleans, optional
            If true, from the low-value to the high-value the values in all
            binned fields will be accumulated.  This is useful for instance
            when adding an unweighted CellMassMsun to a radial plot, as it will
            show mass interior to that radius.  The first value is for the
            x-axis, the second value for the y-axis.  Note that accumulation
            will only be along each row or column.
        x_bins : int, optional
            How many bins should there be in the x-axis variable?
        x_log : boolean, optional
            Should the bin edges be log-spaced?
        x_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        y_bins : int, optional
            How many bins should there be in the y-axis variable?
        y_log : boolean, optional
            Should the bin edges be log-spaced?
        y_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        lazy_reader : boolean, optional
            If this is false, all of the data will be read into memory before
            any processing occurs.  It defaults to true, and grids are binned
            on a one-by-one basis.  Note that parallel computation requires
            this to be true.
        fractional : boolean
            If true, the plot will be normalized to the sum of all the binned
            values.

        Returns
        -------
        plot : `yt.visualization.plot_types.PlotTypes.PhasePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.profiles.BinnedProfile2D : This is the object that does the
                                   transformation of raw data into a 1D
                                   profile.
        
        Examples
        --------
        This will show the mass-distribution in the Density-Temperature plane.
        
        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> reg = pf.h.region([0.1, 0.2, 0.3], [0.0, 0.1, 0.2],
        ...                   [0.2, 0.3, 0.4])
        >>> pc.add_phase_object(reg, ["Density", "Temperature", "CellMassMsun"],
        ...                     weight = None)
        """
        if x_bounds is None:
            x_min, x_max = data_source.quantities["Extrema"](
                                    field_x, non_zero = x_log,
                                    lazy_reader=lazy_reader)[0]
        else:
            x_min, x_max = x_bounds
        if y_bounds is None:
            y_min, y_max = data_source.quantities["Extrema"](
                                    field_y, non_zero = y_log,
                                    lazy_reader=lazy_reader)[0]
        else:
            y_min, y_max = y_bounds
        profile = BinnedProfile2D(data_source,
                                  x_bins, field_x, x_min, x_max, x_log,
                                  y_bins, field_y, y_min, y_max, y_log,
                                  lazy_reader)
        # This is a fallback, in case we forget.
        if field_z.startswith("CellMass") or \
           field_z.startswith("CellVolume"):
            mylog.warning("Setting weight to None")
            weight = None
        profile.add_fields(field_z, weight=weight, accumulation=accumulation, fractional=fractional)
        self._current_field = field_z
        self.profile = profile
        self.scale = {True:'log', False:'linear'}.get(
                data_source.pf.field_info[field_z].take_log, "log")
        self._setup_plot()

    def _setup_plot(self):
        xax = AxisSpec()
        xax.title = self.profile.x_bin_field
        xax.bounds = (self.profile._x_bins[0],
                      self.profile._x_bins[-1])
        xax.scale = {True: 'log', False: 'linear'}[self.profile._x_log]
        xax.calculate_ticks()

        yax = AxisSpec()
        yax.title = self.profile.y_bin_field
        yax.bounds = (self.profile._y_bins[0],
                      self.profile._y_bins[-1])
        yax.scale = {True: 'log', False: 'linear'}[self.profile._y_log]
        yax.calculate_ticks()

        cbar = ColorbarSpec()
        cbar.title = self._current_field
        if self.scale == 'log':
            nz = (self.profile[self._current_field] > 0)
            mi = self.profile[self._current_field][nz].min()
        else:
            mi = self.profile[self._current_field][nz].min()
        ma = self.profile[self._current_field].max()
        cbar.bounds = (mi, ma)
        cbar.cmap = 'algae'
        cbar.scale = self.scale
        cbar.calculate_ticks()

        self.plot = ImagePlotContainer()
        self.plot.image = self.profile[self._current_field].transpose()
        self.plot.x_spec = xax
        self.plot.y_spec = yax
        self.plot.cbar = cbar

class PhasePlotterExtWidget(PhasePlotter):
    _ext_widget_id = None
    _widget_name = "phase_plot"

    def _setup_plot(self):
        if self._ext_widget_id is None: return
        PhasePlotter._setup_plot(self)
        # Now self.plot exists
        from yt.gui.reason.bottle_mods import PayloadHandler
        ph = PayloadHandler()
        # We set up an x axis, y axis, colorbar, and image
        xax = self._convert_axis(self.plot.x_spec)
        yax = self._convert_axis(self.plot.y_spec)
        cbar = self._convert_axis(self.plot.cbar)
        cbar['cmap_image'] = self._get_cbar_image()
        # This is a historical artifact
        raw_data = self.plot.image[::-1,:]

        if self.plot.cbar.scale == 'log':
            func = na.log10
        else:
            func = lambda a: a
        raw_data = na.repeat(raw_data, 3, axis=0)
        raw_data = na.repeat(raw_data, 3, axis=1)
        to_plot = apply_colormap(raw_data, self.plot.cbar.bounds,
                                 self.plot.cbar.cmap, func)
        if self.plot.cbar.scale == 'log':
            # Now we white-out all those regions
            #import pdb;pdb.set_trace()
            to_plot[raw_data == 0.0,:] = 255
        pngs = write_png_to_string(to_plot)
        img_data = base64.b64encode(pngs)
        payload = {'xax':xax, 'yax':yax, 'cbar':cbar,
                   'type': 'widget_payload', 'widget_id': self._ext_widget_id,
                   'image_data': img_data}
        ph.add_payload(payload)

    def _convert_ticks(self, tick_locs, bounds, func, height = 400):
        # height can be a length too; doesn't quite matter.
        mi, ma = func(bounds)
        ticks = []
        for v1,v2 in zip(tick_locs, func(tick_locs)):
            if v2 < mi or v2 > ma: continue
            p = height - height * (v2 - mi)/(ma - mi)
            ticks.append((p,v1,v2))
            #print v1, v2, mi, ma, height, p
        return ticks

    def _convert_axis(self, spec):
        func = lambda a: a
        if spec.scale == 'log': func = na.log10
        tick_info = self._convert_ticks(spec.ticks, spec.bounds, func)
        ax = {'ticks':tick_info,
              'title': spec.title}
        return ax

    def _get_cbar_image(self, height = 400, width = 40):
        # Right now there's just the single 'cmap', but that will eventually
        # change.  I think?
        vals = na.mgrid[1:0:height * 1j] * na.ones(width)[:,None]
        vals = vals.transpose()
        to_plot = apply_colormap(vals)
        pngs = write_png_to_string(to_plot)
        img_data = base64.b64encode(pngs)
        return img_data

class ProfilePlotter(object):
    scale = None
    _current_field = None

    def __init__(self, data_source, field_x, field_y, 
                 weight="CellMassMsun", accumulation=False,
                 x_bins=128, x_log=True, x_bounds=None,
                 lazy_reader=True, fractional=False):
        if x_bounds is None:
            x_min, x_max = data_source.quantities["Extrema"](
                                    field_x, non_zero = x_log,
                                    lazy_reader=lazy_reader)[0]
        else:
            x_min, x_max = x_bounds
        profile = BinnedProfile1D(data_source,
                                  x_bins, field_x, x_min, x_max, x_log,
                                  lazy_reader)
        # This is a fallback, in case we forget.
        if field_y.startswith("CellMass") or \
           field_y.startswith("CellVolume"):
            mylog.warning("Setting weight to None")
            weight = None
        profile.add_fields(field_y, weight=weight, accumulation=accumulation, fractional=fractional)
        self._current_field = field_y
        self.profile = profile
        self.scale = {True:'log', False:'linear'}.get(
                data_source.pf.field_info[field_y].take_log, "log")
        self._setup_plot()

    def _setup_plot(self):
        xax = AxisSpec()
        xax.title = self.profile.bin_field
        xax.bounds = (self.profile._bins[0],
                      self.profile._bins[-1])
        xax.scale = {True: 'log', False: 'linear'}[self.profile._x_log]
        xax.calculate_ticks()

        yax = AxisSpec()
        yax.title = self._current_field
        if self.scale == 'log':
            ny = (self.profile[self._current_field] > 0)
            mi = self.profile[self._current_field][ny].min()
        else:
            mi = self.profile[self._current_field].min()
        ma = self.profile[self._current_field].max()
        yax.bounds = (mi, ma)
        yax.scale = self.scale
        yax.calculate_ticks()

        self.plot = PlotContainer()
        self.plot.x_values = self.profile._bins
        self.plot.y_values = self.profile[self._current_field]
        self.plot.x_spec = xax
        self.plot.y_spec = yax

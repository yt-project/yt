"""
This is a simple mechanism for interfacing with Profile and Phase plots

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
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

from functools import wraps
import numpy as na

from yt.data_objects.profiles import \
    BinnedProfile1D, \
    BinnedProfile2D
from .plot_types import ProfilePlot, PhasePlot
from .loglocator import LogLocator

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
        else:
            raise NotImplementedError
        self.ticks = locator(*self.bounds)

class ColorbarSpec(object):
    title = None
    bounds = None
    scale = None
    cmap = None
    ticks = None

    def calculate_ticks(self):
        if self.scale == 'log':
            locator = LogLocator()
        else:
            raise NotImplementedError
        self.ticks = locator(*self.bounds)


class ImagePlotContainer(object):
    x_spec = None
    y_spec = None
    image = None
    cmap = None
    cbar = None

class PhasePlotter(object):
    scale = None
    _current_field = None

    def __init__(self, data_source, field_x, field_y, field_z,
                 weight="CellMassMsun", accumulation=False,
                 x_bins=128, x_log=True, x_bounds=None,
                 y_bins=128, y_log=True, y_bounds=None,
                 lazy_reader=True, id=None,
                 axes = None, figure = None,
                 fractional=False):
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
            An acceptable colormap.  See either raven.color_maps or
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
        id : int, optional
            If specified, this will be the "semi-unique id" of the resultant
            plot.  This should not be set.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
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
        if field_z == "CellMassMsun": weight = None
        profile.add_fields(field_z, weight=weight, accumulation=accumulation, fractional=fractional)
        self._current_field = field_z
        self.profile = profile
        self.scale = {True:'log', False:'linear'}.get(
                data_source.pf.field_info[field_z], "log")
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
        self.plot.image = self.profile[self._current_field]
        self.plot.x_spec = xax
        self.plot.y_spec = yax
        self.plot.cmap = 'algae'
        self.plot.cbar = cbar

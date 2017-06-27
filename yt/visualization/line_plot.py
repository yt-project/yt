"""
A mechanism for plotting field values along a line through a dataset



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from collections import defaultdict
from yt.funcs import \
    iterable
from yt.units.unit_object import \
    Unit
from yt.units.yt_array import \
    YTArray
from yt.visualization.base_plot_types import \
    PlotMPL
from yt.visualization.plot_container import \
    PlotContainer, \
    PlotDictionary, \
    log_transform, \
    linear_transform, \
    invalidate_plot

class LinePlotDictionary(PlotDictionary):
    def __init__(self, data_source):
        super(LinePlotDictionary, self).__init__(data_source)
        self.known_dimensions = {}

    def _sanitize_dimensions(self, item):
        field = self.data_source._determine_fields(item)[0]
        finfo = self.data_source.ds.field_info[field]
        dimensions = Unit(
            finfo.units, registry=self.data_source.ds.unit_registry).dimensions
        if dimensions not in self.known_dimensions:
            self.known_dimensions[dimensions] = item
            ret_item = item
        else:
            ret_item = self.known_dimensions[dimensions]
        return ret_item

    def __getitem__(self, item):
        ret_item = self._sanitize_dimensions(item)
        return super(LinePlotDictionary, self).__getitem__(ret_item)

    def __setitem__(self, item, value):
        ret_item = self._sanitize_dimensions(item)
        super(LinePlotDictionary, self).__setitem__(ret_item, value)

    def __contains__(self, item):
        ret_item = self._sanitize_dimensions(item)
        return super(LinePlotDictionary, self).__contains__(ret_item)

class LinePlot(PlotContainer):
    r"""
    A class for constructing line plots

    Parameters
    ----------

    ds: :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    fields: string
        The name(s) of the field(s) to be plotted.
    start_point: n-element list, tuple, ndarray, or YTArray
        Contains the coordinates of the first point for constructing the line.
        Must contain n elements where n is the dimensionality of the dataset.
    end_point: n-element list, tuple, ndarray, or YTArray
        Contains the coordinates of the first point for constructing the line.
        Must contain n elements where n is the dimensionality of the dataset.
    figure_size: integer or two-element
    resolution: int
        How many points to sample between start_point and end_point for
        constructing the line plot
    """
    _plot_type = 'line_plot'

    def __init__(self, ds, fields, start_point, end_point, resolution,
                 figure_size=5., fontsize=14., labels=None):
        """
        Sets up figure and axes
        """
        self.start_point = _validate_point(start_point, ds)
        self.end_point = _validate_point(end_point, ds)
        self.resolution = resolution
        self._x_unit = None
        self._y_units = {}
        self._titles = {}

        data_source = ds.all_data()

        self.fields = data_source._determine_fields(fields)
        self.plots = LinePlotDictionary(data_source)
        if labels is None:
            self.labels = {}
        else:
            self.labels = labels

        super(LinePlot, self).__init__(data_source, figure_size, fontsize)

        for f in self.fields:
            if f not in self.labels:
                self.labels[f] = f[1]
            finfo = self.data_source.ds._get_field_info(*f)
            if finfo.take_log:
                self._field_transform[f] = log_transform
            else:
                self._field_transform[f] = linear_transform

        self._setup_plots()

    def add_legend(self, field):
        """Adds a legend to the `LinePlot` instance"""
        plot = self.plots[field]
        plot.axes.legend()

    def _setup_plots(self):
        if self._plot_valid is True:
            return
        dimensions_counter = defaultdict(int)
        for field in self.fields:
            fontscale = self._font_properties._size / 14.
            top_buff_size = 0.3*fontscale

            x_axis_size = 1.2*fontscale
            y_axis_size = 0.9*fontscale
            right_buff_size = 0.2*fontscale

            xbins = np.array([x_axis_size, self.figure_size[0],
                              right_buff_size])
            ybins = np.array([y_axis_size, self.figure_size[1], top_buff_size])

            size = [xbins.sum(), ybins.sum()]

            x_frac_widths = xbins/size[0]
            y_frac_widths = ybins/size[1]

            axrect = (
                x_frac_widths[0],
                y_frac_widths[0],
                x_frac_widths[1],
                y_frac_widths[1],
            )

            try:
                plot = self.plots[field]
            except KeyError:
                plot = PlotMPL(self.figure_size, axrect, None, None)
                self.plots[field] = plot

            plot._set_font_properties(self._font_properties, None)

            x, y = self.ds.coordinates.line_plot(
                field, self.start_point, self.end_point, self.resolution)

            if self._x_unit is None:
                unit_x = x.units
            else:
                unit_x = self._x_unit

            if field in self._y_units:
                unit_y = self._y_units[field]
            else:
                unit_y = y.units

            x = x.to(unit_x)
            y = y.to(unit_y)

            plot.axes.plot(x, y, label=self.labels[field])

            if self._field_transform[field] != linear_transform:
                if (y < 0).any():
                    plot.axes.set_yscale('symlog')
                else:
                    plot.axes.set_yscale('log')

            axes_unit_labels = self._get_axes_unit_labels(unit_x, unit_y)

            finfo = self.ds.field_info[field]

            x_label = r'$\rm{Path\ Length' + axes_unit_labels[0]+'}$'

            finfo = self.ds.field_info[field]
            dimensions = Unit(finfo.units, registry=self.ds.unit_registry).dimensions
            dimensions_counter[dimensions] += 1
            if dimensions_counter[dimensions] > 1:
                y_label = (r'$\rm{Multiple\ Fields}$' + r'$\rm{' +
                           axes_unit_labels[1]+'}$')
            else:
                y_label = (finfo.get_latex_display_name() + r'$\rm{' +
                           axes_unit_labels[1]+'}$')

            plot.axes.set_xlabel(x_label)
            plot.axes.set_ylabel(y_label)

            if field in self._titles:
                plot.axes.set_title(self._titles[field])

    @invalidate_plot
    def set_x_unit(self, unit_name):
        """Set the unit to use along the x-axis

        Parameters
        ----------
        unit_name: str
          The name of the unit to use for the x-axis unit
        """
        self._x_unit = unit_name

    @invalidate_plot
    def set_unit(self, field, unit_name):
        """Set the unit used to plot the field

        Parameters
        ----------
        field: str or field tuple
           The name of the field to set the units for
        unit_name: str
           The name of the unit to use for this field
        """
        self._y_units[self.data_source._determine_fields(field)[0]] = unit_name

    @invalidate_plot
    def annotate_title(self, field, title):
        """Set the unit used to plot the field

        Parameters
        ----------
        field: str or field tuple
           The name of the field to set the units for
        title: str
           The title to use for the plot
        """
        self._titles[self.data_source._determine_fields(field)[0]] = title

def _validate_point(point, ds):
    if not iterable(point):
        raise RuntimeError
    if not isinstance(point, YTArray):
        point = ds.arr(point, 'code_length')
    return point

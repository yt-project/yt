from collections import defaultdict
from typing import Optional

import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm

from yt.funcs import is_sequence, mylog
from yt.units.unit_object import Unit  # type: ignore
from yt.units.yt_array import YTArray
from yt.visualization.plot_container import (
    BaseLinePlot,
    PlotDictionary,
    invalidate_plot,
)


class LineBuffer:
    r"""
    LineBuffer(ds, start_point, end_point, npoints, label = None)

    This takes a data source and implements a protocol for generating a
    'pixelized', fixed-resolution line buffer. In other words, LineBuffer
    takes a starting point, ending point, and number of sampling points and
    can subsequently generate YTArrays of field values along the sample points.

    Parameters
    ----------
    ds : :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object holding the data that can be sampled by the
        LineBuffer
    start_point : n-element list, tuple, ndarray, or YTArray
        Contains the coordinates of the first point for constructing the LineBuffer.
        Must contain n elements where n is the dimensionality of the dataset.
    end_point : n-element list, tuple, ndarray, or YTArray
        Contains the coordinates of the first point for constructing the LineBuffer.
        Must contain n elements where n is the dimensionality of the dataset.
    npoints : int
        How many points to sample between start_point and end_point

    Examples
    --------
    >>> lb = yt.LineBuffer(ds, (0.25, 0, 0), (0.25, 1, 0), 100)
    >>> lb[("all", "u")].max()
    0.11562424257143075 dimensionless

    """

    def __init__(self, ds, start_point, end_point, npoints, label=None):
        self.ds = ds
        self.start_point = _validate_point(start_point, ds, start=True)
        self.end_point = _validate_point(end_point, ds)
        self.npoints = npoints
        self.label = label
        self.data = {}

    def keys(self):
        return self.data.keys()

    def __setitem__(self, item, val):
        self.data[item] = val

    def __getitem__(self, item):
        if item in self.data:
            return self.data[item]
        mylog.info("Making a line buffer with %d points of %s", self.npoints, item)
        self.points, self.data[item] = self.ds.coordinates.pixelize_line(
            item, self.start_point, self.end_point, self.npoints
        )

        return self.data[item]

    def __delitem__(self, item):
        del self.data[item]


class LinePlotDictionary(PlotDictionary):
    def __init__(self, data_source):
        super().__init__(data_source)
        self.known_dimensions = {}

    def _sanitize_dimensions(self, item):
        field = self.data_source._determine_fields(item)[0]
        finfo = self.data_source.ds.field_info[field]
        dimensions = Unit(
            finfo.units, registry=self.data_source.ds.unit_registry
        ).dimensions
        if dimensions not in self.known_dimensions:
            self.known_dimensions[dimensions] = item
        return self.known_dimensions[dimensions]

    def __getitem__(self, item):
        ret_item = self._sanitize_dimensions(item)
        return super().__getitem__(ret_item)

    def __setitem__(self, item, value):
        ret_item = self._sanitize_dimensions(item)
        super().__setitem__(ret_item, value)

    def __contains__(self, item):
        ret_item = self._sanitize_dimensions(item)
        return super().__contains__(ret_item)


class LinePlot(BaseLinePlot):
    r"""
    A class for constructing line plots

    Parameters
    ----------

    ds : :class:`yt.data_objects.static_output.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    fields : string / tuple, or list of strings / tuples
        The name(s) of the field(s) to be plotted.
    start_point : n-element list, tuple, ndarray, or YTArray
        Contains the coordinates of the first point for constructing the line.
        Must contain n elements where n is the dimensionality of the dataset.
    end_point : n-element list, tuple, ndarray, or YTArray
        Contains the coordinates of the first point for constructing the line.
        Must contain n elements where n is the dimensionality of the dataset.
    npoints : int
        How many points to sample between start_point and end_point for
        constructing the line plot
    figure_size : int or two-element iterable of ints
        Size in inches of the image.
        Default: 5 (5x5)
    fontsize : int
        Font size for all text in the plot.
        Default: 14
    field_labels : dictionary
        Keys should be the field names. Values should be latex-formattable
        strings used in the LinePlot legend
        Default: None


    Example
    -------

    >>> import yt

    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> plot = yt.LinePlot(ds, "density", [0, 0, 0], [1, 1, 1], 512)
    >>> plot.add_legend("density")
    >>> plot.set_x_unit("cm")
    >>> plot.set_unit("density", "kg/cm**3")
    >>> plot.save()

    """
    _plot_dict_type = LinePlotDictionary
    _plot_type = "line_plot"

    _default_figure_size = (5.0, 5.0)
    _default_font_size = 14.0

    def __init__(
        self,
        ds,
        fields,
        start_point,
        end_point,
        npoints,
        figure_size=None,
        fontsize: Optional[float] = None,
        field_labels=None,
    ):
        """
        Sets up figure and axes
        """
        line = LineBuffer(ds, start_point, end_point, npoints, label=None)
        self.lines = [line]
        self._initialize_instance(self, ds, fields, figure_size, fontsize, field_labels)
        self._setup_plots()

    @classmethod
    def _initialize_instance(
        cls, obj, ds, fields, figure_size, fontsize, field_labels=None
    ):
        obj._x_unit = None
        obj._titles = {}

        data_source = ds.all_data()

        obj.fields = data_source._determine_fields(fields)
        obj.include_legend = defaultdict(bool)
        super(LinePlot, obj).__init__(
            data_source, figure_size=figure_size, fontsize=fontsize
        )
        if field_labels is None:
            obj.field_labels = {}
        else:
            obj.field_labels = field_labels
        for f in obj.fields:
            if f not in obj.field_labels:
                obj.field_labels[f] = f[1]

    def _get_axrect(self):
        fontscale = self._font_properties._size / self.__class__._default_font_size
        top_buff_size = 0.35 * fontscale

        x_axis_size = 1.35 * fontscale
        y_axis_size = 0.7 * fontscale
        right_buff_size = 0.2 * fontscale

        if is_sequence(self.figure_size):
            figure_size = self.figure_size
        else:
            figure_size = (self.figure_size, self.figure_size)

        xbins = np.array([x_axis_size, figure_size[0], right_buff_size])
        ybins = np.array([y_axis_size, figure_size[1], top_buff_size])

        x_frac_widths = xbins / xbins.sum()
        y_frac_widths = ybins / ybins.sum()

        return (
            x_frac_widths[0],
            y_frac_widths[0],
            x_frac_widths[1],
            y_frac_widths[1],
        )

    @classmethod
    def from_lines(
        cls, ds, fields, lines, figure_size=None, font_size=None, field_labels=None
    ):
        """
        A class method for constructing a line plot from multiple sampling lines

        Parameters
        ----------

        ds : :class:`yt.data_objects.static_output.Dataset`
            This is the dataset object corresponding to the
            simulation output to be plotted.
        fields : field name or list of field names
            The name(s) of the field(s) to be plotted.
        lines : list of :class:`yt.visualization.line_plot.LineBuffer` instances
            The lines from which to sample data
        figure_size : int or two-element iterable of ints
            Size in inches of the image.
            Default: 5 (5x5)
        font_size : int
            Font size for all text in the plot.
            Default: 14
        field_labels : dictionary
            Keys should be the field names. Values should be latex-formattable
            strings used in the LinePlot legend
            Default: None

        Example
        --------
        >>> ds = yt.load(
        ...     "SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e", step=-1
        ... )
        >>> fields = [field for field in ds.field_list if field[0] == "all"]
        >>> lines = [
        ...     yt.LineBuffer(ds, [0.25, 0, 0], [0.25, 1, 0], 100, label="x = 0.25"),
        ...     yt.LineBuffer(ds, [0.5, 0, 0], [0.5, 1, 0], 100, label="x = 0.5"),
        ... ]
        >>> lines.append()

        >>> plot = yt.LinePlot.from_lines(ds, fields, lines)
        >>> plot.save()

        """
        obj = cls.__new__(cls)
        obj.lines = lines
        cls._initialize_instance(obj, ds, fields, figure_size, font_size, field_labels)
        obj._setup_plots()
        return obj

    def _setup_plots(self):
        if self._plot_valid:
            return
        for plot in self.plots.values():
            plot.axes.cla()
        for line in self.lines:
            dimensions_counter = defaultdict(int)
            for field in self.fields:
                finfo = self.ds.field_info[field]
                dimensions = Unit(
                    finfo.units, registry=self.ds.unit_registry
                ).dimensions
                dimensions_counter[dimensions] += 1
            for field in self.fields:
                # get plot instance
                plot = self._get_plot_instance(field)

                # calculate x and y
                x, y = self.ds.coordinates.pixelize_line(
                    field, line.start_point, line.end_point, line.npoints
                )

                # scale x and y to proper units
                if self._x_unit is None:
                    unit_x = x.units
                else:
                    unit_x = self._x_unit

                unit_y = plot.norm_handler.display_units

                x.convert_to_units(unit_x)
                y.convert_to_units(unit_y)

                # determine legend label
                str_seq = []
                str_seq.append(line.label)
                str_seq.append(self.field_labels[field])
                delim = "; "
                legend_label = delim.join(filter(None, str_seq))

                # apply plot to matplotlib axes
                plot.axes.plot(x, y, label=legend_label)

                # apply log transforms if requested
                norm = plot.norm_handler.get_norm(data=y)
                y_norm_type = type(norm)
                if y_norm_type is Normalize:
                    plot.axes.set_yscale("linear")
                elif y_norm_type is LogNorm:
                    plot.axes.set_yscale("log")
                elif y_norm_type is SymLogNorm:
                    plot.axes.set_yscale("symlog")
                else:
                    raise NotImplementedError(
                        f"LinePlot doesn't support y norm with type {type(norm)}"
                    )

                # set font properties
                plot._set_font_properties(self._font_properties, None)

                # set x and y axis labels
                axes_unit_labels = self._get_axes_unit_labels(unit_x, unit_y)

                if self._xlabel is not None:
                    x_label = self._xlabel
                else:
                    x_label = r"$\rm{Path\ Length" + axes_unit_labels[0] + "}$"

                if self._ylabel is not None:
                    y_label = self._ylabel
                else:
                    finfo = self.ds.field_info[field]
                    dimensions = Unit(
                        finfo.units, registry=self.ds.unit_registry
                    ).dimensions
                    if dimensions_counter[dimensions] > 1:
                        y_label = (
                            r"$\rm{Multiple\ Fields}$"
                            + r"$\rm{"
                            + axes_unit_labels[1]
                            + "}$"
                        )
                    else:
                        y_label = (
                            finfo.get_latex_display_name()
                            + r"$\rm{"
                            + axes_unit_labels[1]
                            + "}$"
                        )

                plot.axes.set_xlabel(x_label)
                plot.axes.set_ylabel(y_label)

                # apply title
                if field in self._titles:
                    plot.axes.set_title(self._titles[field])

                # apply legend
                dim_field = self.plots._sanitize_dimensions(field)
                if self.include_legend[dim_field]:
                    plot.axes.legend()

        self._plot_valid = True

    @invalidate_plot
    def annotate_legend(self, field):
        """
        Adds a legend to the `LinePlot` instance. The `_sanitize_dimensions`
        call ensures that a legend label will be added for every field of
        a multi-field plot
        """
        dim_field = self.plots._sanitize_dimensions(field)
        self.include_legend[dim_field] = True

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
    def set_unit(self, field, new_unit):
        """Set the unit used to plot the field

        Parameters
        ----------
        field: str or field tuple
           The name of the field to set the units for
        new_unit: string or Unit object
        """
        field = self.data_source._determine_fields(field)[0]
        pnh = self.plots[field].norm_handler
        pnh.display_units = new_unit

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


def _validate_point(point, ds, start=False):
    if not is_sequence(point):
        raise RuntimeError("Input point must be array-like")
    if not isinstance(point, YTArray):
        point = ds.arr(point, "code_length", dtype=np.float64)
    if len(point.shape) != 1:
        raise RuntimeError("Input point must be a 1D array")
    if point.shape[0] < ds.dimensionality:
        raise RuntimeError("Input point must have an element for each dimension")
    # need to pad to 3D elements to avoid issues later
    if point.shape[0] < 3:
        if start:
            val = 0
        else:
            val = 1
        point = np.append(point.d, [val] * (3 - ds.dimensionality)) * point.uq
    return point

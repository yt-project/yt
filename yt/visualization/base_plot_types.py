import sys
import warnings
from abc import ABC
from io import BytesIO
from typing import TYPE_CHECKING, Optional, Tuple, Union

import matplotlib
import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.ticker import LogFormatterMathtext
from packaging.version import Version

from yt.funcs import get_interactivity, is_sequence, matplotlib_style_context, mylog
from yt.visualization._handlers import ColorbarHandler, NormHandler

from ._commons import (
    MPL_VERSION,
    get_canvas,
    get_symlog_majorticks,
    get_symlog_minorticks,
    validate_image_name,
)

if TYPE_CHECKING:
    from matplotlib.axis import Axis
    from matplotlib.figure import Figure

BACKEND_SPECS = {
    "GTK": ["backend_gtk", "FigureCanvasGTK", "FigureManagerGTK"],
    "GTKAgg": ["backend_gtkagg", "FigureCanvasGTKAgg", None],
    "GTKCairo": ["backend_gtkcairo", "FigureCanvasGTKCairo", None],
    "MacOSX": ["backend_macosx", "FigureCanvasMac", "FigureManagerMac"],
    "Qt4Agg": ["backend_qt4agg", "FigureCanvasQTAgg", None],
    "Qt5Agg": ["backend_qt5agg", "FigureCanvasQTAgg", None],
    "TkAgg": ["backend_tkagg", "FigureCanvasTkAgg", None],
    "WX": ["backend_wx", "FigureCanvasWx", None],
    "WXAgg": ["backend_wxagg", "FigureCanvasWxAgg", None],
    "GTK3Cairo": [
        "backend_gtk3cairo",
        "FigureCanvasGTK3Cairo",
        "FigureManagerGTK3Cairo",
    ],
    "GTK3Agg": ["backend_gtk3agg", "FigureCanvasGTK3Agg", "FigureManagerGTK3Agg"],
    "WebAgg": ["backend_webagg", "FigureCanvasWebAgg", None],
    "nbAgg": ["backend_nbagg", "FigureCanvasNbAgg", "FigureManagerNbAgg"],
    "agg": ["backend_agg", "FigureCanvasAgg", None],
}


class CallbackWrapper:
    def __init__(self, viewer, window_plot, frb, field, font_properties, font_color):
        self.frb = frb
        self.data = frb.data_source
        self._axes = window_plot.axes
        self._figure = window_plot.figure
        if len(self._axes.images) > 0:
            self.raw_image_shape = self._axes.images[0]._A.shape
            if viewer._has_swapped_axes:
                # store the original un-transposed shape
                self.raw_image_shape = self.raw_image_shape[1], self.raw_image_shape[0]
        if frb.axis < 3:
            DD = frb.ds.domain_width
            xax = frb.ds.coordinates.x_axis[frb.axis]
            yax = frb.ds.coordinates.y_axis[frb.axis]
            self._period = (DD[xax], DD[yax])
        self.ds = frb.ds
        self.xlim = viewer.xlim
        self.ylim = viewer.ylim
        self._swap_axes = viewer._has_swapped_axes
        self._flip_horizontal = viewer._flip_horizontal  # needed for quiver
        self._flip_vertical = viewer._flip_vertical  # needed for quiver
        # an important note on _swap_axes: _swap_axes will swap x,y arguments
        # in callbacks (e.g., plt.plot(x,y) will be plt.plot(y, x). The xlim
        # and ylim arguments above, and internal callback references to coordinates
        # are the **unswapped** ranges.
        self._axes_unit_names = viewer._axes_unit_names
        if "OffAxisSlice" in viewer._plot_type:
            self._type_name = "CuttingPlane"
        else:
            self._type_name = viewer._plot_type
        self.aspect = window_plot._aspect
        self.font_properties = font_properties
        self.font_color = font_color
        self.field = field


class PlotMPL:
    """A base class for all yt plots made using matplotlib, that is backend independent."""

    def __init__(
        self,
        fsize,
        axrect,
        *,
        norm_handler: NormHandler,
        figure: Optional["Figure"] = None,
        axes: Optional["Axis"] = None,
    ):
        """Initialize PlotMPL class"""
        import matplotlib.figure

        self._plot_valid = True
        if figure is None:
            if not is_sequence(fsize):
                fsize = (fsize, fsize)
            self.figure = matplotlib.figure.Figure(figsize=fsize, frameon=True)
        else:
            figure.set_size_inches(fsize)
            self.figure = figure
        if axes is None:
            self._create_axes(axrect)
        else:
            axes.cla()
            axes.set_position(axrect)
            self.axes = axes
        self.interactivity = get_interactivity()

        figure_canvas, figure_manager = self._get_canvas_classes()
        self.canvas = figure_canvas(self.figure)
        if figure_manager is not None:
            self.manager = figure_manager(self.canvas, 1)

        self.axes.tick_params(
            which="both", axis="both", direction="in", top=True, right=True
        )

        self.norm_handler = norm_handler

    def _create_axes(self, axrect):
        self.axes = self.figure.add_axes(axrect)

    def _get_canvas_classes(self):

        if self.interactivity:
            key = str(matplotlib.get_backend())
        else:
            key = "agg"

        try:
            module, fig_canvas, fig_manager = BACKEND_SPECS[key]
        except KeyError:
            return

        mod = __import__(
            "matplotlib.backends",
            globals(),
            locals(),
            [module],
            0,
        )
        submod = getattr(mod, module)
        FigureCanvas = getattr(submod, fig_canvas)
        if fig_manager is not None:
            FigureManager = getattr(submod, fig_manager)
            return FigureCanvas, FigureManager

        return FigureCanvas, None

    def save(self, name, mpl_kwargs=None, canvas=None):
        """Choose backend and save image to disk"""

        if mpl_kwargs is None:
            mpl_kwargs = {}
        if "papertype" not in mpl_kwargs and MPL_VERSION < Version("3.3.0"):
            mpl_kwargs["papertype"] = "auto"

        name = validate_image_name(name)

        try:
            canvas = get_canvas(self.figure, name)
        except ValueError:
            canvas = self.canvas

        mylog.info("Saving plot %s", name)
        with matplotlib_style_context():
            canvas.print_figure(name, **mpl_kwargs)
        return name

    def show(self):
        try:
            self.manager.show()
        except AttributeError:
            self.canvas.show()

    def _get_labels(self):
        ax = self.axes
        labels = ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels()
        labels += ax.xaxis.get_minorticklabels()
        labels += ax.yaxis.get_minorticklabels()
        labels += [
            ax.title,
            ax.xaxis.label,
            ax.yaxis.label,
            ax.xaxis.get_offset_text(),
            ax.yaxis.get_offset_text(),
        ]
        return labels

    def _set_font_properties(self, font_properties, font_color):
        for label in self._get_labels():
            label.set_fontproperties(font_properties)
            if font_color is not None:
                label.set_color(font_color)

    def _repr_png_(self):
        from ._mpl_imports import FigureCanvasAgg

        canvas = FigureCanvasAgg(self.figure)
        f = BytesIO()
        with matplotlib_style_context():
            canvas.print_figure(f)
        f.seek(0)
        return f.read()


class ImagePlotMPL(PlotMPL, ABC):
    """A base class for yt plots made using imshow"""

    _default_font_size = 18.0

    def __init__(
        self,
        fsize=None,
        axrect=None,
        caxrect=None,
        *,
        norm_handler: NormHandler,
        colorbar_handler: ColorbarHandler,
        figure: Optional["Figure"] = None,
        axes: Optional["Axis"] = None,
        cax: Optional["Axis"] = None,
    ):
        """Initialize ImagePlotMPL class object"""
        self.colorbar_handler = colorbar_handler
        _missing_layout_specs = [_ is None for _ in (fsize, axrect, caxrect)]

        if all(_missing_layout_specs):
            fsize, axrect, caxrect = self._get_best_layout()
        elif any(_missing_layout_specs):
            raise TypeError(
                "ImagePlotMPL cannot be initialized with partially specified layout."
            )

        super().__init__(
            fsize, axrect, norm_handler=norm_handler, figure=figure, axes=axes
        )

        if cax is None:
            self.cax = self.figure.add_axes(caxrect)
        else:
            cax.cla()
            cax.set_position(caxrect)
            self.cax = cax

    def _setup_layout_constraints(
        self, figure_size: Union[Tuple[float, float], float], fontsize: float
    ):
        # Setup base layout attributes
        # derived classes need to call this before super().__init__
        # but they are free to do other stuff in between

        if isinstance(figure_size, tuple):
            assert len(figure_size) == 2
            assert all(isinstance(_, float) for _ in figure_size)
            self._figure_size = figure_size
        else:
            assert isinstance(figure_size, float)
            self._figure_size = (figure_size, figure_size)

        self._draw_axes = True
        fontscale = float(fontsize) / self.__class__._default_font_size
        if fontscale < 1.0:
            fontscale = np.sqrt(fontscale)

        self._cb_size = 0.0375 * self._figure_size[0]
        self._ax_text_size = [1.2 * fontscale, 0.9 * fontscale]
        self._top_buff_size = 0.30 * fontscale
        self._aspect = 1.0

    def _reset_layout(self) -> None:
        size, axrect, caxrect = self._get_best_layout()
        self.axes.set_position(axrect)
        self.cax.set_position(caxrect)
        self.figure.set_size_inches(*size)

    def _init_image(self, data, extent, aspect):
        """Store output of imshow in image variable"""

        if MPL_VERSION < Version("3.2"):
            # with MPL 3.1 we use np.inf as a mask instead of np.nan
            # this is done in CoordinateHandler.sanitize_buffer_fill_values
            # however masking with inf is problematic here when we search for the max value
            # so here we revert to nan
            # see https://github.com/yt-project/yt/pull/2517 and https://github.com/yt-project/yt/pull/3793
            data[~np.isfinite(data)] = np.nan

        norm = self.norm_handler.get_norm(data)
        extent = [float(e) for e in extent]

        if self._transform is None:
            # sets the transform to be an ax.TransData object, where the
            # coordinate system of the data is controlled by the xlim and ylim
            # of the data.
            transform = self.axes.transData
        else:
            transform = self._transform

        self._validate_axes_extent(extent, transform)

        self.image = self.axes.imshow(
            data.to_ndarray(),
            origin="lower",
            extent=extent,
            norm=norm,
            aspect=aspect,
            cmap=self.colorbar_handler.cmap,
            interpolation="nearest",
            transform=transform,
        )
        self._set_axes(norm)

    def _set_axes(self, norm: Normalize) -> None:
        if isinstance(norm, SymLogNorm):
            formatter = LogFormatterMathtext(linthresh=norm.linthresh)
            self.cb = self.figure.colorbar(self.image, self.cax, format=formatter)
            self.cb.set_ticks(
                get_symlog_majorticks(
                    linthresh=norm.linthresh, vmin=norm.vmin, vmax=norm.vmax
                )
            )
        else:
            self.cb = self.figure.colorbar(self.image, self.cax)
        self.cax.tick_params(which="both", axis="y", direction="in")

        fmt_kwargs = dict(style="scientific", scilimits=(-2, 3), useMathText=True)
        self.image.axes.ticklabel_format(**fmt_kwargs)
        if type(norm) not in (LogNorm, SymLogNorm):
            try:
                self.cb.ax.ticklabel_format(**fmt_kwargs)
            except AttributeError as exc:
                if MPL_VERSION < Version("3.5.0"):
                    warnings.warn(
                        "Failed to format colorbar ticks. "
                        "This is expected when using the set_norm method "
                        "with some matplotlib classes (e.g. TwoSlopeNorm) "
                        "with matplotlib versions older than 3.5\n"
                        "Please try upgrading matplotlib to a more recent version. "
                        "If the problem persists, please file a report to "
                        "https://github.com/yt-project/yt/issues/new"
                    )
                else:
                    raise exc
        if self.colorbar_handler.draw_minorticks:
            if isinstance(norm, SymLogNorm):
                if Version("3.2.0") <= MPL_VERSION < Version("3.5.0b"):
                    # no known working method to draw symlog minor ticks
                    # see https://github.com/yt-project/yt/issues/3535
                    pass
                else:
                    flinthresh = 10 ** np.floor(np.log10(norm.linthresh))
                    absmax = np.abs((norm.vmin, norm.vmax)).max()
                    if (absmax - flinthresh) / absmax < 0.1:
                        flinthresh /= 10
                    mticks = get_symlog_minorticks(flinthresh, norm.vmin, norm.vmax)
                    if MPL_VERSION < Version("3.5.0b"):
                        # https://github.com/matplotlib/matplotlib/issues/21258
                        mticks = self.image.norm(mticks)
                    self.cax.yaxis.set_ticks(mticks, minor=True)

            elif isinstance(norm, LogNorm):
                self.cax.minorticks_on()
                self.cax.xaxis.set_visible(False)

            else:
                self.cax.minorticks_on()
        else:
            self.cax.minorticks_off()

        self.image.axes.set_facecolor(self.colorbar_handler.background_color)

    def _validate_axes_extent(self, extent, transform):
        # if the axes are cartopy GeoAxes, this checks that the axes extent
        # is properly set.

        if "cartopy" not in sys.modules:
            # cartopy isn't already loaded, nothing to do here
            return

        from cartopy.mpl.geoaxes import GeoAxes

        if isinstance(self.axes, GeoAxes):
            # some projections have trouble when passing extents at or near the
            # limits. So we only set_extent when the plot is a subset of the
            # globe, within the tolerance of the transform.

            # note that `set_extent` here is setting the extent of the axes.
            # still need to pass the extent arg to imshow in order to
            # ensure that it is properly scaled. also note that set_extent
            # expects values in the coordinates of the transform: it will
            # calculate the coordinates in the projection.
            global_extent = transform.x_limits + transform.y_limits
            thresh = transform.threshold
            if all(
                abs(extent[ie]) < (abs(global_extent[ie]) - thresh) for ie in range(4)
            ):
                self.axes.set_extent(extent, crs=transform)

    def _get_best_layout(self):
        # this method is called in ImagePlotMPL.__init__
        # required attributes
        # - self._figure_size: Union[float, Tuple[float, float]]
        # - self._aspect: float
        # - self._ax_text_size: Tuple[float, float]
        # - self._draw_axes: bool
        # - self.colorbar_handler: ColorbarHandler

        # optional attribtues
        # - self._unit_aspect: float

        # Ensure the figure size along the long axis is always equal to _figure_size
        unit_aspect = getattr(self, "_unit_aspect", 1)
        if is_sequence(self._figure_size):
            x_fig_size, y_fig_size = self._figure_size
            y_fig_size *= unit_aspect
        else:
            x_fig_size = y_fig_size = self._figure_size
            scaling = self._aspect / unit_aspect
            if scaling < 1:
                x_fig_size *= scaling
            else:
                y_fig_size /= scaling

        if self.colorbar_handler.draw_cbar:
            cb_size = self._cb_size
            cb_text_size = self._ax_text_size[1] + 0.45
        else:
            cb_size = x_fig_size * 0.04
            cb_text_size = 0.0

        if self._draw_axes:
            x_axis_size = self._ax_text_size[0]
            y_axis_size = self._ax_text_size[1]
        else:
            x_axis_size = x_fig_size * 0.04
            y_axis_size = y_fig_size * 0.04

        top_buff_size = self._top_buff_size

        if not self._draw_axes and not self.colorbar_handler.draw_cbar:
            x_axis_size = 0.0
            y_axis_size = 0.0
            cb_size = 0.0
            cb_text_size = 0.0
            top_buff_size = 0.0

        xbins = np.array([x_axis_size, x_fig_size, cb_size, cb_text_size])
        ybins = np.array([y_axis_size, y_fig_size, top_buff_size])

        size = [xbins.sum(), ybins.sum()]

        x_frac_widths = xbins / size[0]
        y_frac_widths = ybins / size[1]

        # axrect is the rectangle defining the area of the
        # axis object of the plot.  Its range goes from 0 to 1 in
        # x and y directions.  The first two values are the x,y
        # start values of the axis object (lower left corner), and the
        # second two values are the size of the axis object.  To get
        # the upper right corner, add the first x,y to the second x,y.
        axrect = (
            x_frac_widths[0],
            y_frac_widths[0],
            x_frac_widths[1],
            y_frac_widths[1],
        )

        # caxrect is the rectangle defining the area of the colorbar
        # axis object of the plot.  It is defined just as the axrect
        # tuple is.
        caxrect = (
            x_frac_widths[0] + x_frac_widths[1],
            y_frac_widths[0],
            x_frac_widths[2],
            y_frac_widths[1],
        )

        return size, axrect, caxrect

    def _toggle_axes(self, choice, draw_frame=None):
        """
        Turn on/off displaying the axis ticks and labels for a plot.

        Parameters
        ----------
        choice : boolean
            If True, set the axes to be drawn. If False, set the axes to not be
            drawn.
        """
        self._draw_axes = choice
        self._draw_frame = draw_frame
        if draw_frame is None:
            draw_frame = choice
        if self.colorbar_handler.has_background_color and not draw_frame:
            # workaround matplotlib's behaviour
            # last checked with Matplotlib 3.5
            warnings.warn(
                f"Previously set background color {self.colorbar_handler.background_color} "
                "has no effect. Pass `draw_frame=True` if you wish to preserve background color.",
                stacklevel=4,
            )
        self.axes.set_frame_on(draw_frame)
        self.axes.get_xaxis().set_visible(choice)
        self.axes.get_yaxis().set_visible(choice)
        self._reset_layout()

    def _toggle_colorbar(self, choice: bool):
        """
        Turn on/off displaying the colorbar for a plot

        choice = True or False
        """
        self.colorbar_handler.draw_cbar = choice
        self.cax.set_visible(choice)
        size, axrect, caxrect = self._get_best_layout()
        self.axes.set_position(axrect)
        self.cax.set_position(caxrect)
        self.figure.set_size_inches(*size)

    def _get_labels(self):
        labels = super()._get_labels()
        cbax = self.cb.ax
        labels += cbax.yaxis.get_ticklabels()
        labels += [cbax.yaxis.label, cbax.yaxis.get_offset_text()]
        return labels

    def hide_axes(self, *, draw_frame=None):
        """
        Hide the axes for a plot including ticks and labels
        """
        self._toggle_axes(False, draw_frame)
        return self

    def show_axes(self):
        """
        Show the axes for a plot including ticks and labels
        """
        self._toggle_axes(True)
        return self

    def hide_colorbar(self):
        """
        Hide the colorbar for a plot including ticks and labels
        """
        self._toggle_colorbar(False)
        return self

    def show_colorbar(self):
        """
        Show the colorbar for a plot including ticks and labels
        """
        self._toggle_colorbar(True)
        return self


def get_multi_plot(nx, ny, colorbar="vertical", bw=4, dpi=300, cbar_padding=0.4):
    r"""Construct a multiple axes plot object, with or without a colorbar, into
    which multiple plots may be inserted.

    This will create a set of :class:`matplotlib.axes.Axes`, all lined up into
    a grid, which are then returned to the user and which can be used to plot
    multiple plots on a single figure.

    Parameters
    ----------
    nx : int
        Number of axes to create along the x-direction
    ny : int
        Number of axes to create along the y-direction
    colorbar : {'vertical', 'horizontal', None}, optional
        Should Axes objects for colorbars be allocated, and if so, should they
        correspond to the horizontal or vertical set of axes?
    bw : number
        The base height/width of an axes object inside the figure, in inches
    dpi : number
        The dots per inch fed into the Figure instantiation

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        The figure created inside which the axes reside
    tr : list of list of :class:`matplotlib.axes.Axes` objects
        This is a list, where the inner list is along the x-axis and the outer
        is along the y-axis
    cbars : list of :class:`matplotlib.axes.Axes` objects
        Each of these is an axes onto which a colorbar can be placed.

    Notes
    -----
    This is a simple implementation for a common use case.  Viewing the source
    can be instructive, and is encouraged to see how to generate more
    complicated or more specific sets of multiplots for your own purposes.
    """
    import matplotlib.figure

    hf, wf = 1.0 / ny, 1.0 / nx
    fudge_x = fudge_y = 1.0
    if colorbar is None:
        fudge_x = fudge_y = 1.0
    elif colorbar.lower() == "vertical":
        fudge_x = nx / (cbar_padding + nx)
        fudge_y = 1.0
    elif colorbar.lower() == "horizontal":
        fudge_x = 1.0
        fudge_y = ny / (cbar_padding + ny)
    fig = matplotlib.figure.Figure((bw * nx / fudge_x, bw * ny / fudge_y), dpi=dpi)
    from ._mpl_imports import FigureCanvasAgg

    fig.set_canvas(FigureCanvasAgg(fig))
    fig.subplots_adjust(
        wspace=0.0, hspace=0.0, top=1.0, bottom=0.0, left=0.0, right=1.0
    )
    tr = []
    for j in range(ny):
        tr.append([])
        for i in range(nx):
            left = i * wf * fudge_x
            bottom = fudge_y * (1.0 - (j + 1) * hf) + (1.0 - fudge_y)
            ax = fig.add_axes([left, bottom, wf * fudge_x, hf * fudge_y])
            tr[-1].append(ax)
    cbars = []
    if colorbar is None:
        pass
    elif colorbar.lower() == "horizontal":
        for i in range(nx):
            # left, bottom, width, height
            # Here we want 0.10 on each side of the colorbar
            # We want it to be 0.05 tall
            # And we want a buffer of 0.15
            ax = fig.add_axes(
                [
                    wf * (i + 0.10) * fudge_x,
                    hf * fudge_y * 0.20,
                    wf * (1 - 0.20) * fudge_x,
                    hf * fudge_y * 0.05,
                ]
            )
            cbars.append(ax)
    elif colorbar.lower() == "vertical":
        for j in range(ny):
            ax = fig.add_axes(
                [
                    wf * (nx + 0.05) * fudge_x,
                    hf * fudge_y * (ny - (j + 0.95)),
                    wf * fudge_x * 0.05,
                    hf * fudge_y * 0.90,
                ]
            )
            ax.clear()
            cbars.append(ax)
    return fig, tr, cbars

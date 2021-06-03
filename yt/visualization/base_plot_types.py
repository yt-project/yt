from distutils.version import LooseVersion
from io import BytesIO

import matplotlib
import numpy as np

from yt.funcs import (
    get_brewer_cmap,
    get_interactivity,
    is_sequence,
    matplotlib_style_context,
    mylog,
)

from ._commons import get_canvas, validate_image_name

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
            self.image = self._axes.images[0]
        if frb.axis < 3:
            DD = frb.ds.domain_width
            xax = frb.ds.coordinates.x_axis[frb.axis]
            yax = frb.ds.coordinates.y_axis[frb.axis]
            self._period = (DD[xax], DD[yax])
        self.ds = frb.ds
        self.xlim = viewer.xlim
        self.ylim = viewer.ylim
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

    def __init__(self, fsize, axrect, figure, axes):
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

        for which in ["major", "minor"]:
            for axis in "xy":
                self.axes.tick_params(
                    which=which, axis=axis, direction="in", top=True, right=True
                )

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
        if "papertype" not in mpl_kwargs and LooseVersion(
            matplotlib.__version__
        ) < LooseVersion("3.3.0"):
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
                label.set_color(self.font_color)

    def _repr_png_(self):
        from ._mpl_imports import FigureCanvasAgg

        canvas = FigureCanvasAgg(self.figure)
        f = BytesIO()
        with matplotlib_style_context():
            canvas.print_figure(f)
        f.seek(0)
        return f.read()


class ImagePlotMPL(PlotMPL):
    """A base class for yt plots made using imshow"""

    def __init__(self, fsize, axrect, caxrect, zlim, figure, axes, cax):
        """Initialize ImagePlotMPL class object"""
        super().__init__(fsize, axrect, figure, axes)
        self.zmin, self.zmax = zlim
        if cax is None:
            self.cax = self.figure.add_axes(caxrect)
        else:
            cax.cla()
            cax.set_position(caxrect)
            self.cax = cax

    def _init_image(self, data, cbnorm, cblinthresh, cmap, extent, aspect):
        """Store output of imshow in image variable"""
        cbnorm_kwargs = dict(
            vmin=float(self.zmin) if self.zmin is not None else None,
            vmax=float(self.zmax) if self.zmax is not None else None,
        )
        if cbnorm == "log10":
            cbnorm_cls = matplotlib.colors.LogNorm
        elif cbnorm == "linear":
            cbnorm_cls = matplotlib.colors.Normalize
        elif cbnorm == "symlog":
            # if cblinthresh is not specified, try to come up with a reasonable default
            vmin = float(np.nanmin(data))
            vmax = float(np.nanmax(data))
            if cblinthresh is None:
                cblinthresh = np.nanmin(np.absolute(data)[data != 0])

            cbnorm_kwargs.update(dict(linthresh=cblinthresh, vmin=vmin, vmax=vmax))
            MPL_VERSION = LooseVersion(matplotlib.__version__)
            if MPL_VERSION >= "3.2.0":
                # note that this creates an inconsistency between mpl versions
                # since the default value previous to mpl 3.4.0 is np.e
                # but it is only exposed since 3.2.0
                cbnorm_kwargs["base"] = 10

            cbnorm_cls = matplotlib.colors.SymLogNorm
        else:
            raise ValueError(f"Unknown value `cbnorm` == {cbnorm}")

        norm = cbnorm_cls(**cbnorm_kwargs)

        extent = [float(e) for e in extent]
        # tuple colormaps are from palettable (or brewer2mpl)
        if isinstance(cmap, tuple):
            cmap = get_brewer_cmap(cmap)

        if self._transform is None:
            # sets the transform to be an ax.TransData object, where the
            # coordiante system of the data is controlled by the xlim and ylim
            # of the data.
            transform = self.axes.transData
        else:
            transform = self._transform
        if hasattr(self.axes, "set_extent"):
            # CartoPy hangs if we do not set_extent before imshow if we are
            # displaying a small subset of the globe.  What I believe happens is
            # that the transform for the points on the outside results in
            # infinities, and then the scipy.spatial cKDTree hangs trying to
            # identify nearest points.
            #
            # Also, set_extent is defined by cartopy, so not all axes will have
            # it as a method.
            #
            # A potential downside is that other images may change, but I believe
            # the result of imshow is to set_extent *regardless*.  This just
            # changes the order in which it happens.
            #
            # NOTE: This is currently commented out because it breaks in some
            # instances.  It is left as a historical note because we will
            # eventually need some form of it.
            # self.axes.set_extent(extent)
            pass
        self.image = self.axes.imshow(
            data.to_ndarray(),
            origin="lower",
            extent=extent,
            norm=norm,
            aspect=aspect,
            cmap=cmap,
            interpolation="nearest",
            transform=transform,
        )
        if cbnorm == "symlog":
            formatter = matplotlib.ticker.LogFormatterMathtext(linthresh=cblinthresh)
            self.cb = self.figure.colorbar(self.image, self.cax, format=formatter)
            if np.nanmin(data) >= 0.0:
                yticks = [np.nanmin(data).v] + list(
                    10
                    ** np.arange(
                        np.rint(np.log10(cblinthresh)),
                        np.ceil(np.log10(np.nanmax(data))) + 1,
                    )
                )
            elif np.nanmax(data) <= 0.0:
                yticks = (
                    list(
                        -(
                            10
                            ** np.arange(
                                np.floor(np.log10(-np.nanmin(data))),
                                np.rint(np.log10(cblinthresh)) - 1,
                                -1,
                            )
                        )
                    )
                    + [np.nanmax(data).v]
                )
            else:
                yticks = (
                    list(
                        -(
                            10
                            ** np.arange(
                                np.floor(np.log10(-np.nanmin(data))),
                                np.rint(np.log10(cblinthresh)) - 1,
                                -1,
                            )
                        )
                    )
                    + [0]
                    + list(
                        10
                        ** np.arange(
                            np.rint(np.log10(cblinthresh)),
                            np.ceil(np.log10(np.nanmax(data))) + 1,
                        )
                    )
                )
            self.cb.set_ticks(yticks)
        else:
            self.cb = self.figure.colorbar(self.image, self.cax)
        for which in ["major", "minor"]:
            self.cax.tick_params(which=which, axis="y", direction="in")

    def _get_best_layout(self):

        # Ensure the figure size along the long axis is always equal to _figure_size
        if is_sequence(self._figure_size):
            x_fig_size = self._figure_size[0]
            y_fig_size = self._figure_size[1]
        else:
            x_fig_size = self._figure_size
            y_fig_size = self._figure_size / self._aspect

        if hasattr(self, "_unit_aspect"):
            y_fig_size = y_fig_size * self._unit_aspect

        if self._draw_colorbar:
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

        if not self._draw_axes and not self._draw_colorbar:
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
        if draw_frame is None:
            draw_frame = choice
        self._draw_axes = choice
        self._draw_frame = draw_frame
        self.axes.set_frame_on(draw_frame)
        self.axes.get_xaxis().set_visible(choice)
        self.axes.get_yaxis().set_visible(choice)
        size, axrect, caxrect = self._get_best_layout()
        self.axes.set_position(axrect)
        self.cax.set_position(caxrect)
        self.figure.set_size_inches(*size)

    def _toggle_colorbar(self, choice):
        """
        Turn on/off displaying the colorbar for a plot

        choice = True or False
        """
        self._draw_colorbar = choice
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

    def hide_axes(self, draw_frame=None):
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

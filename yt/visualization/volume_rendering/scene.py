import builtins
import functools
from collections import OrderedDict
from typing import List, Optional

import numpy as np

from yt.config import ytcfg
from yt.funcs import mylog
from yt.units.dimensions import length  # type: ignore
from yt.units.unit_registry import UnitRegistry  # type: ignore
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.exceptions import YTNotInsideNotebook
from yt.visualization._commons import get_canvas, validate_image_name

from .camera import Camera
from .render_source import (
    BoxSource,
    CoordinateVectorSource,
    GridSource,
    LineSource,
    MeshSource,
    OpaqueSource,
    PointSource,
    RenderSource,
    VolumeSource,
)
from .zbuffer_array import ZBuffer


class Scene:

    """A virtual landscape for a volume rendering.

    The Scene class is meant to be the primary container for the
    new volume rendering framework. A single scene may contain
    several Camera and RenderSource instances, and is the primary
    driver behind creating a volume rendering.

    This sets up the basics needed to add sources and cameras.
    This does very little setup, and requires additional input
    to do anything useful.

    Examples
    --------

    This example shows how to create an empty scene and add a VolumeSource
    and a Camera.

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import (
    ...     Camera, Scene, create_volume_source)
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> sc = Scene()
    >>> source = create_volume_source(ds.all_data(), "density")
    >>> sc.add_source(source)
    >>> cam = sc.add_camera()
    >>> im = sc.render()

    Alternatively, you can use the create_scene function to set up defaults
    and then modify the Scene later:

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> sc = yt.create_scene(ds)
    >>> # Modify camera, sources, etc...
    >>> im = sc.render()

    """

    _current = None
    _camera = None
    _unit_registry = None

    def __init__(self):
        r"""Create a new Scene instance"""
        super().__init__()
        self.sources = OrderedDict()
        self._last_render = None
        # A non-public attribute used to get around the fact that we can't
        # pass kwargs into _repr_png_()
        self._sigma_clip = None

    def get_source(self, source_num=0):
        """Returns the volume rendering source indexed by ``source_num``"""
        return list(self.sources.values())[source_num]

    def __getitem__(self, item):
        if item in self.sources:
            return self.sources[item]
        return self.get_source(item)

    @property
    def opaque_sources(self):
        """
        Iterate over opaque RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.items():
            if isinstance(source, OpaqueSource) or issubclass(
                OpaqueSource, type(source)
            ):
                yield k, source

    @property
    def transparent_sources(self):
        """
        Iterate over transparent RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.items():
            if not isinstance(source, OpaqueSource):
                yield k, source

    def add_source(self, render_source, keyname=None):
        """Add a render source to the scene.

        This will autodetect the type of source.

        Parameters
        ----------
        render_source:
            :class:`yt.visualization.volume_rendering.render_source.RenderSource`
            A source to contribute to the volume rendering scene.

        keyname: string (optional)
            The dictionary key used to reference the source in the sources
            dictionary.
        """
        if keyname is None:
            keyname = "source_%02i" % len(self.sources)

        data_sources = (VolumeSource, MeshSource, GridSource)

        if isinstance(render_source, data_sources):
            self._set_new_unit_registry(render_source.data_source.ds.unit_registry)

        line_annotation_sources = (GridSource, BoxSource, CoordinateVectorSource)

        if isinstance(render_source, line_annotation_sources):
            lens_str = str(self.camera.lens)
            if "fisheye" in lens_str or "spherical" in lens_str:
                raise NotImplementedError(
                    "Line annotation sources are not supported for %s."
                    % (type(self.camera.lens).__name__),
                )

        if isinstance(render_source, (LineSource, PointSource)):
            if isinstance(render_source.positions, YTArray):
                render_source.positions = (
                    self.arr(render_source.positions).in_units("code_length").d
                )

        self.sources[keyname] = render_source

        return self

    def __setitem__(self, key, value):
        return self.add_source(value, key)

    def _set_new_unit_registry(self, input_registry):
        self.unit_registry = UnitRegistry(
            add_default_symbols=False, lut=input_registry.lut
        )

        # Validate that the new unit registry makes sense
        current_scaling = self.unit_registry["unitary"][0]
        if current_scaling != input_registry["unitary"][0]:
            for source in self.sources.items():
                data_source = getattr(source, "data_source", None)
                if data_source is None:
                    continue
                scaling = data_source.ds.unit_registry["unitary"][0]
                if scaling != current_scaling:
                    raise NotImplementedError(
                        "Simultaneously rendering data from datasets with "
                        "different units is not supported"
                    )

    def render(self, camera=None):
        r"""Render all sources in the Scene.

        Use the current state of the Scene object to render all sources
        currently in the scene.  Returns the image array.  If you want to
        save the output to a file, call the save() function.

        Parameters
        ----------
        camera: :class:`Camera`, optional
            If specified, use a different :class:`Camera` to render the scene.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` instance containing
        the current rendering image.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> im = sc.render()
        >>> sc.save(sigma_clip=4.0, render=False)

        Altneratively, if you do not need the image array, you can just call
        ``save`` as follows.

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> sc.save(sigma_clip=4.0)

        """
        mylog.info("Rendering scene (Can take a while).")
        if camera is None:
            camera = self.camera
        assert camera is not None
        self._validate()
        bmp = self.composite(camera=camera)
        self._last_render = bmp
        return bmp

    def _render_on_demand(self, render):
        # checks for existing render before rendering, in most cases we want to
        # render every time, but in some cases pulling the previous render is
        # desirable (e.g., if only changing sigma_clip or
        # saving after a call to sc.show()).

        if self._last_render is not None and not render:
            mylog.info("Found previously rendered image to save.")
            return

        if self._last_render is None:
            mylog.warning("No previously rendered image found, rendering now.")
        elif render:
            mylog.warning(
                "Previously rendered image exists, but rendering anyway. "
                "Supply 'render=False' to save previously rendered image directly."
            )
        self.render()

    def _get_render_sources(self):
        return [s for s in self.sources.values() if isinstance(s, RenderSource)]

    def _setup_save(self, fname, render) -> str:

        self._render_on_demand(render)

        rensources = self._get_render_sources()
        if fname is None:
            # if a volume source present, use its affiliated ds for fname
            if len(rensources) > 0:
                rs = rensources[0]
                basename = rs.data_source.ds.basename
                if isinstance(rs.field, str):
                    field = rs.field
                else:
                    field = rs.field[-1]
                fname = f"{basename}_Render_{field}"
            # if no volume source present, use a default filename
            else:
                fname = "Render_opaque"

        fname = validate_image_name(fname)
        mylog.info("Saving rendered image to %s", fname)
        return fname

    def save(
        self,
        fname: Optional[str] = None,
        sigma_clip: Optional[float] = None,
        render: bool = True,
    ):
        r"""Saves a rendered image of the Scene to disk.

        Once you have created a scene, this saves an image array to disk with
        an optional filename. This function calls render() to generate an
        image array, unless the render parameter is set to False, in which case
        the most recently rendered scene is used if it exists.

        Parameters
        ----------
        fname: string, optional
            If specified, save the rendering as to the file "fname".
            If unspecified, it creates a default based on the dataset filename.
            The file format is inferred from the filename's suffix.
            Supported formats depend on which version of matplotlib is installed.

            Default: None
        sigma_clip: float, optional
            Image values greater than this number times the standard deviation
            plus the mean of the image will be clipped before saving. Useful
            for enhancing images as it gets rid of rare high pixel values.
            Default: None

            floor(vals > std_dev*sigma_clip + mean)
        render: boolean, optional
            If True, will always render the scene before saving.
            If False, will use results of previous render if it exists.
            Default: True

        Returns
        -------
            Nothing

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> sc.save("test.png", sigma_clip=4)

        When saving multiple images without modifying the scene (camera,
        sources,etc.), render=False can be used to avoid re-rendering.
        This is useful for generating images at a range of sigma_clip values:

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> # save with different sigma clipping values
        >>> sc.save("raw.png")  # The initial render call happens here
        >>> sc.save("clipped_2.png", sigma_clip=2, render=False)
        >>> sc.save("clipped_4.png", sigma_clip=4, render=False)

        """
        fname = self._setup_save(fname, render)

        # We can render pngs natively but for other formats we defer to
        # matplotlib.
        if fname.endswith(".png"):
            self._last_render.write_png(fname, sigma_clip=sigma_clip)
        else:
            from matplotlib.figure import Figure

            shape = self._last_render.shape
            fig = Figure((shape[0] / 100.0, shape[1] / 100.0))
            canvas = get_canvas(fig, fname)

            ax = fig.add_axes([0, 0, 1, 1])
            ax.set_axis_off()
            out = self._last_render
            if sigma_clip is not None:
                max_val = out._clipping_value(sigma_clip)
            else:
                max_val = out[:, :, :3].max()
            alpha = 255 * out[:, :, 3].astype("uint8")
            out = np.clip(out[:, :, :3] / max_val, 0.0, 1.0) * 255
            out = np.concatenate([out.astype("uint8"), alpha[..., None]], axis=-1)
            # not sure why we need rot90, but this makes the orientation
            # match the png writer
            ax.imshow(np.rot90(out), origin="lower")
            canvas.print_figure(fname, dpi=100)

    def save_annotated(
        self,
        fname: Optional[str] = None,
        label_fmt: Optional[str] = None,
        text_annotate=None,
        dpi: int = 100,
        sigma_clip: Optional[float] = None,
        render: bool = True,
        tf_rect: Optional[List[float]] = None,
    ):
        r"""Saves the most recently rendered image of the Scene to disk,
        including an image of the transfer function and and user-defined
        text.

        Once you have created a scene and rendered that scene to an image
        array, this saves that image array to disk with an optional filename.
        If an image has not yet been rendered for the current scene object,
        it forces one and writes it out.

        Parameters
        ----------
        fname: string, optional
            If specified, save the rendering as a bitmap to the file "fname".
            If unspecified, it creates a default based on the dataset filename.
            Default: None
        sigma_clip: float, optional
            Image values greater than this number times the standard deviation
            plus the mean of the image will be clipped before saving. Useful
            for enhancing images as it gets rid of rare high pixel values.
            Default: None

            floor(vals > std_dev*sigma_clip + mean)
        dpi: integer, optional
            By default, the resulting image will be the same size as the camera
            parameters.  If you supply a dpi, then the image will be scaled
            accordingly (from the default 100 dpi)
        label_fmt : str, optional
           A format specifier (e.g., label_fmt="%.2g") to use in formatting
           the data values that label the transfer function colorbar.
        text_annotate : list of iterables
           Any text that you wish to display on the image.  This should be an
           list containing a tuple of coordinates (in normalized figure
           coordinates), the text to display, and, optionally, a dictionary of
           keyword/value pairs to pass through to the matplotlib text()
           function.

           Each item in the main list is a separate string to write.
        render: boolean, optional
            If True, will render the scene before saving.
            If False, will use results of previous render if it exists.
            Default: True
        tf_rect : sequence of floats, optional
           A rectangle that defines the location of the transfer
           function legend.  This is only used for the case where
           there are multiple volume sources with associated transfer
           functions.  tf_rect is of the form [x0, y0, width, height],
           in figure coordinates.

        Returns
        -------
            Nothing


        Examples
        --------

        >>> sc.save_annotated(
        ...     "fig.png",
        ...     text_annotate=[
        ...         [
        ...             (0.05, 0.05),
        ...             f"t = {ds.current_time.d}",
        ...             dict(horizontalalignment="left"),
        ...         ],
        ...         [
        ...             (0.5, 0.95),
        ...             "simulation title",
        ...             dict(color="y", fontsize="24", horizontalalignment="center"),
        ...         ],
        ...     ],
        ... )

        """
        fname = self._setup_save(fname, render)

        ax = self._show_mpl(
            self._last_render.swapaxes(0, 1), sigma_clip=sigma_clip, dpi=dpi
        )

        # number of transfer functions?
        num_trans_func = 0
        for rs in self._get_render_sources():
            if hasattr(rs, "transfer_function"):
                num_trans_func += 1

        # which transfer function?
        if num_trans_func == 1:
            rs = self._get_render_sources()[0]
            tf = rs.transfer_function
            label = rs.data_source.ds._get_field_info(rs.field).get_label()
            self._annotate(ax.axes, tf, rs, label=label, label_fmt=label_fmt)
        else:
            # set the origin and width and height of the colorbar region
            if tf_rect is None:
                tf_rect = [0.80, 0.12, 0.12, 0.9]
            cbx0, cby0, cbw, cbh = tf_rect

            cbh_each = cbh / num_trans_func

            for i, rs in enumerate(self._get_render_sources()):
                ax = self._render_figure.add_axes(
                    [cbx0, cby0 + i * cbh_each, 0.8 * cbw, 0.8 * cbh_each]
                )
                try:
                    tf = rs.transfer_function
                except AttributeError:
                    pass
                else:
                    label = rs.data_source.ds._get_field_info(rs.field).get_label()
                    self._annotate_multi(ax, tf, rs, label=label, label_fmt=label_fmt)

        # any text?
        if text_annotate is not None:
            f = self._render_figure
            for t in text_annotate:
                xy = t[0]
                string = t[1]
                if len(t) == 3:
                    opt = t[2]
                else:
                    opt = dict()

                # sane default
                if "color" not in opt:
                    opt["color"] = "w"

                ax.axes.text(xy[0], xy[1], string, transform=f.transFigure, **opt)

        self._render_figure.canvas = get_canvas(self._render_figure, fname)
        self._render_figure.tight_layout()
        self._render_figure.savefig(fname, facecolor="black", pad_inches=0)

    def _show_mpl(self, im, sigma_clip=None, dpi=100):
        from matplotlib.figure import Figure

        s = im.shape
        self._render_figure = Figure(figsize=(s[1] / float(dpi), s[0] / float(dpi)))
        self._render_figure.clf()
        ax = self._render_figure.add_subplot(111)
        ax.set_position([0, 0, 1, 1])

        if sigma_clip is not None:
            nim = im / im._clipping_value(sigma_clip)
            nim[nim > 1.0] = 1.0
            nim[nim < 0.0] = 0.0
        else:
            nim = im
        axim = ax.imshow(nim[:, :, :3] / nim[:, :, :3].max(), interpolation="bilinear")

        return axim

    def _annotate(self, ax, tf, source, label="", label_fmt=None):
        ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_visible(False)
        ax.get_yaxis().set_ticks([])
        cb = self._render_figure.colorbar(
            ax.images[0], pad=0.0, fraction=0.05, drawedges=True
        )
        tf.vert_cbar(
            ax=cb.ax,
            label=label,
            label_fmt=label_fmt,
            resolution=self.camera.resolution[0],
            log_scale=source.log_field,
        )

    def _annotate_multi(self, ax, tf, source, label, label_fmt):
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        tf.vert_cbar(
            ax=ax,
            label=label,
            label_fmt=label_fmt,
            resolution=self.camera.resolution[0],
            log_scale=source.log_field,
            size=6,
        )

    def _validate(self):
        r"""Validate the current state of the scene."""

        for source in self.sources.values():
            source._validate()
        return

    def composite(self, camera=None):
        r"""Create a composite image of the current scene.

        First iterate over the opaque sources and set the ZBuffer.
        Then iterate over the transparent sources, rendering from the value
        of the zbuffer to the front of the box. Typically this function is
        accessed through the .render() command.

        Parameters
        ----------
        camera: :class:`Camera`, optional
            If specified, use a specific :class:`Camera` to render the scene.

        Returns
        -------
        im: :class:`ImageArray`
            ImageArray instance of the current rendering image.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> im = sc.composite()

        """
        if camera is None:
            camera = self.camera
        empty = camera.lens.new_image(camera)
        opaque = ZBuffer(empty, np.full(empty.shape[:2], np.inf))

        for _, source in self.opaque_sources:
            source.render(camera, zbuffer=opaque)
            im = source.zbuffer.rgba

        for _, source in self.transparent_sources:
            im = source.render(camera, zbuffer=opaque)
            opaque.rgba = im

        # rotate image 180 degrees so orientation agrees with e.g.
        # a PlotWindow plot
        return np.rot90(im, k=2)

    def add_camera(self, data_source=None, lens_type="plane-parallel", auto=False):
        r"""Add a new camera to the Scene.

        The camera is defined by a position (the location of the camera
        in the simulation domain,), a focus (the point at which the
        camera is pointed), a width (the width of the snapshot that will
        be taken, a resolution (the number of pixels in the image), and
        a north_vector (the "up" direction in the resulting image). A
        camera can use a variety of different Lens objects.

        If the scene already has a camera associated with it, this function
        will create a new camera and discard the old one.

        Parameters
        ----------
        data_source: :class:`AMR3DData` or :class:`Dataset`, optional
            This is the source to be rendered, which can be any arbitrary yt
            data object or dataset.
        lens_type: string, optional
            This specifies the type of lens to use for rendering. Current
            options are 'plane-parallel', 'perspective', and 'fisheye'. See
            :class:`yt.visualization.volume_rendering.lens.Lens` for details.
            Default: 'plane-parallel'
        auto: boolean
            If True, build smart defaults using the data source extent. This
            can be time-consuming to iterate over the entire dataset to find
            the positional bounds. Default: False

        Examples
        --------

        In this example, the camera is set using defaults that are chosen
        to be reasonable for the argument Dataset.

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Camera, Scene
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> sc = Scene()
        >>> sc.add_camera()

        Here, we set the camera properties manually:

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Camera, Scene
        >>> sc = Scene()
        >>> cam = sc.add_camera()
        >>> cam.position = np.array([0.5, 0.5, -1.0])
        >>> cam.focus = np.array([0.5, 0.5, 0.0])
        >>> cam.north_vector = np.array([1.0, 0.0, 0.0])

        Finally, we create a camera with a non-default lens:

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> sc = Scene()
        >>> sc.add_camera(ds, lens_type="perspective")

        """
        self._camera = Camera(self, data_source, lens_type, auto)
        return self.camera

    @property
    def camera(self):
        r"""The camera property.

        This is the default camera that will be used when rendering. Can be set
        manually, but Camera type will be checked for validity.
        """
        return self._camera

    @camera.setter
    def camera(self, value):
        value.width = self.arr(value.width)
        value.focus = self.arr(value.focus)
        value.position = self.arr(value.position)
        self._camera = value

    @camera.deleter
    def camera(self):
        del self._camera
        self._camera = None

    @property
    def unit_registry(self):
        ur = self._unit_registry
        if ur is None:
            ur = UnitRegistry()
            # This will be updated when we add a volume source
            ur.add("unitary", 1.0, length)
        self._unit_registry = ur
        return self._unit_registry

    @unit_registry.setter
    def unit_registry(self, value):
        self._unit_registry = value
        if self.camera is not None:
            self.camera.width = YTArray(
                self.camera.width.in_units("unitary"), registry=value
            )
            self.camera.focus = YTArray(
                self.camera.focus.in_units("unitary"), registry=value
            )
            self.camera.position = YTArray(
                self.camera.position.in_units("unitary"), registry=value
            )

    @unit_registry.deleter
    def unit_registry(self):
        del self._unit_registry
        self._unit_registry = None

    def set_camera(self, camera):
        r"""

        Set the camera to be used by this scene.

        """
        self.camera = camera

    def get_camera(self):
        r"""

        Get the camera currently used by this scene.

        """
        return self.camera

    def annotate_domain(self, ds, color=None):
        r"""

        Modifies this scene by drawing the edges of the computational domain.
        This adds a new BoxSource to the scene corresponding to the domain
        boundaries and returns the modified scene object.

        Parameters
        ----------

        ds : :class:`yt.data_objects.static_output.Dataset`
            This is the dataset object corresponding to the
            simulation being rendered. Used to get the domain bounds.
        color : array_like of shape (4,), optional
            The RGBA value to use to draw the domain boundaries.
            Default is black with an alpha of 1.0.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> sc.annotate_domain(ds)
        >>> im = sc.render()

        """
        box_source = BoxSource(ds.domain_left_edge, ds.domain_right_edge, color=color)
        self.add_source(box_source)
        return self

    def annotate_grids(
        self, data_source, alpha=0.3, cmap=None, min_level=None, max_level=None
    ):
        r"""

        Modifies this scene by drawing the edges of the AMR grids.
        This adds a new GridSource to the scene that represents the AMR grid
        and returns the resulting Scene object.

        Parameters
        ----------

        data_source: :class:`~yt.data_objects.api.DataContainer`
            The data container that will be used to identify grids to draw.
        alpha : float
            The opacity of the grids to draw.
        cmap : color map name
            The color map to use to map resolution levels to color.
        min_level : int, optional
            Minimum level to draw
        max_level : int, optional
            Maximum level to draw


        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> sc.annotate_grids(ds.all_data())
        >>> im = sc.render()

        """
        if cmap is None:
            cmap = ytcfg.get("yt", "default_colormap")
        grids = GridSource(
            data_source,
            alpha=alpha,
            cmap=cmap,
            min_level=min_level,
            max_level=max_level,
        )
        self.add_source(grids)
        return self

    def annotate_mesh_lines(self, color=None, alpha=1.0):
        """

        Modifies this Scene by drawing the mesh line boundaries
        on all MeshSources.

        Parameters
        ----------
        color : array_like of shape (4,), optional
            The RGBA value to use to draw the mesh lines.
            Default is black with an alpha of 1.0.
        alpha : float, optional
            The opacity of the mesh lines. Default is 255 (solid).

        """
        for _, source in self.opaque_sources:
            if isinstance(source, MeshSource):
                source.annotate_mesh_lines(color=color, alpha=alpha)
        return self

    def annotate_axes(self, colors=None, alpha=1.0):
        r"""

        Modifies this scene by drawing the coordinate axes.
        This adds a new CoordinateVectorSource to the scene
        and returns the modified scene object.

        Parameters
        ----------
        colors: array-like of shape (3,4), optional
            The RGBA values to use to draw the x, y, and z vectors. The default
            is  [[1, 0, 0, alpha], [0, 1, 0, alpha], [0, 0, 1, alpha]] where
            ``alpha`` is set by the parameter below. If ``colors`` is set then
            ``alpha`` is ignored.
        alpha : float, optional
            The opacity of the vectors.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> sc.annotate_axes(alpha=0.5)
        >>> im = sc.render()

        """
        coords = CoordinateVectorSource(colors, alpha)
        self.add_source(coords)
        return self

    def show(self, sigma_clip=None):
        r"""This will send the most recently rendered image to the IPython
        notebook.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will send the current image of this Scene
        to the notebook for display. If there is no current image, it will
        run the render() method on this Scene before sending the result to the
        notebook.

        If yt can't determine if it's inside an IPython session, this will raise
        YTNotInsideNotebook.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

        >>> sc = yt.create_scene(ds)
        >>> sc.show()

        """
        if "__IPYTHON__" in dir(builtins):
            from IPython.display import display

            self._sigma_clip = sigma_clip
            display(self)
        else:
            raise YTNotInsideNotebook

    _arr = None

    @property
    def arr(self):
        """Converts an array into a :class:`yt.units.yt_array.YTArray`

        The returned YTArray will be dimensionless by default, but can be
        cast to arbitrary units using the ``units`` keyword argument.

        Parameters
        ----------

        input_array : Iterable
            A tuple, list, or array to attach units to
        units: String unit specification, unit symbol object, or astropy units object
            The units of the array. Powers must be specified using python syntax
            (cm**3, not cm^3).
        dtype : string or NumPy dtype object
            The dtype of the returned array data

        Examples
        --------

        >>> a = sc.arr([1, 2, 3], "cm")
        >>> b = sc.arr([4, 5, 6], "m")
        >>> a + b
        YTArray([ 401.,  502.,  603.]) cm
        >>> b + a
        YTArray([ 4.01,  5.02,  6.03]) m

        Arrays returned by this function know about the scene's unit system

        >>> a = sc.arr(np.ones(5), "unitary")
        >>> a.in_units("Mpc")
        YTArray([ 1.00010449,  1.00010449,  1.00010449,  1.00010449,
                 1.00010449]) Mpc

        """
        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(YTArray, registry=self.unit_registry)
        return self._arr

    _quan = None

    @property
    def quan(self):
        """Converts an scalar into a :class:`yt.units.yt_array.YTQuantity`

        The returned YTQuantity will be dimensionless by default, but can be
        cast to arbitrary units using the ``units`` keyword argument.

        Parameters
        ----------

        input_scalar : an integer or floating point scalar
            The scalar to attach units to
        units : String unit specification, unit symbol object, or astropy
            units
        input_units : deprecated in favor of 'units'
            The units of the quantity. Powers must be specified using python
            syntax (cm**3, not cm^3).
        dtype : string or NumPy dtype object
            The dtype of the array data.

        Examples
        --------

        >>> a = sc.quan(1, "cm")
        >>> b = sc.quan(2, "m")
        >>> a + b
        201.0 cm
        >>> b + a
        2.01 m

        Quantities created this way automatically know about the unit system
        of the scene

        >>> a = ds.quan(5, "unitary")
        >>> a.in_cgs()
        1.543e+25 cm

        """
        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(YTQuantity, registry=self.unit_registry)
        return self._quan

    def _repr_png_(self):
        if self._last_render is None:
            self.render()
        png = self._last_render.write_png(
            filename=None, sigma_clip=self._sigma_clip, background="black"
        )
        self._sigma_clip = None
        return png

    def __repr__(self):
        disp = "<Scene Object>:"
        disp += "\nSources: \n"
        for k, v in self.sources.items():
            disp += f"    {k}: {v}\n"
        disp += "Camera: \n"
        disp += f"    {self.camera}"
        return disp

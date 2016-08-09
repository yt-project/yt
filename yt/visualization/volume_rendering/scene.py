"""
The volume rendering Scene class.

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


import functools
import numpy as np
from collections import OrderedDict
from yt.config import \
    ytcfg
from yt.funcs import mylog, get_image_suffix
from yt.extern.six import iteritems, itervalues, string_types
from yt.units.dimensions import \
    length
from yt.units.unit_registry import \
    UnitRegistry
from yt.units.yt_array import \
    YTQuantity, \
    YTArray
from .camera import Camera
from .render_source import \
    OpaqueSource, \
    BoxSource, \
    CoordinateVectorSource, \
    GridSource, \
    RenderSource, \
    MeshSource, \
    VolumeSource, \
    PointSource, \
    LineSource
from .zbuffer_array import ZBuffer
from yt.extern.six.moves import builtins
from yt.utilities.exceptions import YTNotInsideNotebook

class Scene(object):

    """A virtual landscape for a volume rendering.

    The Scene class is meant to be the primary container for the
    new volume rendering framework. A single scene may contain
    several Camera and RenderSource instances, and is the primary
    driver behind creating a volume rendering.

    This sets up the basics needed to add sources and cameras.
    This does very little setup, and requires additional input
    to do anything useful.

    Parameters
    ----------
    None

    Examples
    --------

    This example shows how to create an empty scene and add a VolumeSource
    and a Camera.

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import Scene, VolumeSource, Camera
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> sc = Scene()
    >>> source = VolumeSource(ds.all_data(), 'density')
    >>> sc.add_source(source)
    >>> cam = sc.add_camera()
    >>> im = sc.render()

    Alternatively, you can use the create_scene function to set up defaults
    and then modify the Scene later:

    >>> import yt
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>>
    >>> sc = yt.create_scene(ds)
    >>> # Modify camera, sources, etc...
    >>> im = sc.render()

    """

    _current = None
    _camera = None
    _unit_registry = None

    def __init__(self):
        r"""Create a new Scene instance"""
        super(Scene, self).__init__()
        self.sources = OrderedDict()
        self._last_render = None
        # A non-public attribute used to get around the fact that we can't
        # pass kwargs into _repr_png_()
        self._sigma_clip = None

    def get_source(self, source_num):
        """Returns the volume rendering source indexed by ``source_num``"""
        return list(itervalues(self.sources))[source_num]

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
        for k, source in iteritems(self.sources):
            if isinstance(source, OpaqueSource) or \
                    issubclass(OpaqueSource, type(source)):
                yield k, source

    @property
    def transparent_sources(self):
        """
        Iterate over transparent RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in iteritems(self.sources):
            if not isinstance(source, OpaqueSource):
                yield k, source

    def add_source(self, render_source, keyname=None):
        """Add a render source to the scene.

        This will autodetect the type of source.

        Parameters
        ----------
        render_source: an instance of :class:`yt.visualization.volume_rendering.render_source.RenderSource`
            A source to contribute to the volume rendering scene.

        keyname: string (optional)
            The dictionary key used to reference the source in the sources
            dictionary.
        """
        if keyname is None:
            keyname = 'source_%02i' % len(self.sources)

        data_sources = (VolumeSource, MeshSource, GridSource)

        if isinstance(render_source, data_sources):
            self._set_new_unit_registry(
                render_source.data_source.ds.unit_registry)

        line_annotation_sources = (GridSource, BoxSource, CoordinateVectorSource)

        if isinstance(render_source, line_annotation_sources):
            lens_str = str(self.camera.lens)
            if 'fisheye' in lens_str or 'spherical' in lens_str:
                raise NotImplementedError(
                    "Line annotation sources are not supported for %s."
                    % (type(self.camera.lens).__name__), )

        if isinstance(render_source, (LineSource, PointSource)):
            if isinstance(render_source.positions, YTArray):
                render_source.positions = \
                    self.arr(render_source.positions).in_units('code_length').d

        self.sources[keyname] = render_source

        return self

    def __setitem__(self, key, value):
        return self.add_source(value, key)

    def _set_new_unit_registry(self, input_registry):
        self.unit_registry = UnitRegistry(
            add_default_symbols=False,
            lut=input_registry.lut)

        # Validate that the new unit registry makes sense
        current_scaling = self.unit_registry['unitary'][0]
        if current_scaling != input_registry['unitary'][0]:
            for source in self.sources.items():
                data_source = getattr(source, 'data_source', None)
                if data_source is None:
                    continue
                scaling = data_source.ds.unit_registry['unitary'][0]
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
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> im = sc.render()
        >>> sc.save(sigma_clip=4.0)

        """
        mylog.info("Rendering scene (Can take a while).")
        if camera is None:
            camera = self.camera
        assert(camera is not None)
        self._validate()
        bmp = self.composite(camera=camera)
        self._last_render = bmp
        return bmp

    def save(self, fname=None, sigma_clip=None):
        r"""Saves the most recently rendered image of the Scene to disk.

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

        Returns
        -------
            Nothing

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> sc.render()
        >>> sc.save('test.png', sigma_clip=4)

        Or alternatively:

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
        >>> sc = yt.create_scene(ds)
        >>> # save with different sigma clipping values
        >>> sc.save('raw.png')
        >>> sc.save('clipped_2.png', sigma_clip=2)
        >>> sc.save('clipped_4.png', sigma_clip=4)

        """
        if fname is None:
            sources = list(itervalues(self.sources))
            rensources = [s for s in sources if isinstance(s, RenderSource)]
            # if a volume source present, use its affiliated ds for fname
            if len(rensources) > 0:
                rs = rensources[0]
                basename = rs.data_source.ds.basename
                if isinstance(rs.field, string_types):
                    field = rs.field
                else:
                    field = rs.field[-1]
                fname = "%s_Render_%s.png" % (basename, field)
            # if no volume source present, use a default filename
            else:
                fname = "Render_opaque.png"
        suffix = get_image_suffix(fname)
        if suffix == '':
            suffix = '.png'
            fname = '%s%s' % (fname, suffix)

        self.render()

        mylog.info("Saving render %s", fname)
        self._last_render.write_png(fname, sigma_clip=sigma_clip)


    def save_annotated(self, fname=None, label_fmt=None,
                       text_annotate=None, dpi=100, sigma_clip=None):
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


        Returns
        -------
            Nothing


        Examples
        --------

        >>> sc.save_annotated("fig.png", 
        ...                   text_annotate=[[(0.05, 0.05),
        ...                                   "t = {}".format(ds.current_time.d),
        ...                                   dict(horizontalalignment="left")],
        ...                                  [(0.5,0.95),
        ...                                   "simulation title",
        ...                                   dict(color="y", fontsize="24",
        ...                                        horizontalalignment="center")]])

        """
        import matplotlib.pyplot as plt

        sources = list(itervalues(self.sources))
        rensources = [s for s in sources if isinstance(s, RenderSource)]

        if fname is None:
            # if a volume source present, use its affiliated ds for fname
            if len(rensources) > 0:
                rs = rensources[0]
                basename = rs.data_source.ds.basename
                if isinstance(rs.field, string_types):
                    field = rs.field
                else:
                    field = rs.field[-1]
                fname = "%s_Render_%s.png" % (basename, field)
            # if no volume source present, use a default filename
            else:
                fname = "Render_opaque.png"
        suffix = get_image_suffix(fname)
        if suffix == '':
            suffix = '.png'
            fname = '%s%s' % (fname, suffix)

        self.render()

        # which transfer function?
        rs = rensources[0]
        tf = rs.transfer_function
        label = rs.data_source.ds._get_field_info(rs.field).get_label()
        if rs.log_field:
            label = r'$\rm{log}\ $' + label

        ax = self._show_mpl(self._last_render.swapaxes(0, 1),
                            sigma_clip=sigma_clip, dpi=dpi)
        self._annotate(ax.axes, tf, rs, label=label, label_fmt=label_fmt)
        plt.tight_layout()

        # any text?
        if text_annotate is not None:
            f = plt.gcf()
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

                plt.text(xy[0], xy[1], string,
                         transform=f.transFigure, **opt)

        plt.savefig(fname, facecolor='black', pad_inches=0)

    def _show_mpl(self, im, sigma_clip=None, dpi=100):
        import matplotlib.pyplot as plt
        s = im.shape
        self._render_figure = plt.figure(1, figsize=(s[1]/float(dpi), s[0]/float(dpi)))
        self._render_figure.clf()
        ax = plt.gca()
        ax.set_position([0, 0, 1, 1])

        if sigma_clip is not None:
            nz = im[im > 0.0]
            nim = im / (nz.mean() + sigma_clip * np.std(nz))
            nim[nim > 1.0] = 1.0
            nim[nim < 0.0] = 0.0
            del nz
        else:
            nim = im
        axim = plt.imshow(nim[:,:,:3]/nim[:,:,:3].max(),
                          interpolation="bilinear")

        return axim

    def _annotate(self, ax, tf, source, label="", label_fmt=None):
        import matplotlib.pyplot as plt
        ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_visible(False)
        ax.get_yaxis().set_ticks([])
        cb = plt.colorbar(ax.images[0], pad=0.0, fraction=0.05,
                          drawedges=True, shrink=0.75)
        tf.vert_cbar(ax=cb.ax, label=label, label_fmt=label_fmt,
                     resolution=self.camera.resolution[0],
                     log_scale=source.log_field)

    def _validate(self):
        r"""Validate the current state of the scene."""

        for k, source in iteritems(self.sources):
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
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
        >>> sc = yt.create_scene(ds)
        >>> # Modify camera, sources, etc...
        >>> im = sc.composite()

        """
        if camera is None:
            camera = self.camera
        empty = camera.lens.new_image(camera)
        opaque = ZBuffer(empty, np.full(empty.shape[:2], np.inf))

        for k, source in self.opaque_sources:
            source.render(camera, zbuffer=opaque)
            im = source.zbuffer.rgba

        for k, source in self.transparent_sources:
            im = source.render(camera, zbuffer=opaque)

        return im

    def add_camera(self, data_source=None, lens_type='plane-parallel',
                   auto=False):
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
        >>> from yt.visualization.volume_rendering.api import Scene, Camera
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> sc = Scene()
        >>> sc.add_camera()

        Here, we set the camera properties manually:

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Scene, Camera
        >>> sc = Scene()
        >>> cam = sc.add_camera()
        >>> cam.position = np.array([0.5, 0.5, -1.0])
        >>> cam.focus = np.array([0.5, 0.5, 0.0])
        >>> cam.north_vector = np.array([1.0, 0.0, 0.0])

        Finally, we create a camera with a non-default lens:

        >>> import yt
        >>> from yt.visualization.volume_rendering.api import Camera
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> sc = Scene()
        >>> sc.add_camera(ds, lens_type='perspective')

        """
        self._camera = Camera(self, data_source, lens_type, auto)
        return self.camera

    def camera():
        doc = r"""The camera property.

        This is the default camera that will be used when rendering. Can be set
        manually, but Camera type will be checked for validity.
        """

        def fget(self):
            return self._camera

        def fset(self, value):
            value.width = self.arr(value.width)
            value.focus = self.arr(value.focus)
            value.position = self.arr(value.position)
            self._camera = value

        def fdel(self):
            del self._camera
            self._camera = None
        return locals()
    camera = property(**camera())

    def unit_registry():
        def fget(self):
            ur = self._unit_registry
            if ur is None:
                ur = UnitRegistry()
                # This will be updated when we add a volume source
                ur.add("unitary", 1.0, length)
            self._unit_registry = ur
            return self._unit_registry

        def fset(self, value):
            self._unit_registry = value
            if self.camera is not None:
                self.camera.width = YTArray(
                    self.camera.width.in_units('unitary'), registry=value)
                self.camera.focus = YTArray(
                    self.camera.focus.in_units('unitary'), registry=value)
                self.camera.position = YTArray(
                    self.camera.position.in_units('unitary'), registry=value)

        def fdel(self):
            del self._unit_registry
            self._unit_registry = None
        return locals()
    unit_registry = property(**unit_registry())

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

        ds : :class:`yt.data_objects.api.Dataset`
            This is the dataset object corresponding to the
            simulation being rendered. Used to get the domain bounds.


        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
        >>> sc = yt.create_scene(ds)
        >>> sc.annotate_domain(ds)
        >>> im = sc.render()

        """
        box_source = BoxSource(ds.domain_left_edge,
                               ds.domain_right_edge,
                               color=color)
        self.add_source(box_source)
        return self

    def annotate_grids(self, data_source, alpha=0.3, cmap=None,
                       min_level=None, max_level=None):
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
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
        >>> sc = yt.create_scene(ds)
        >>> sc.annotate_grids(ds.all_data())
        >>> im = sc.render()

        """
        if cmap is None:
            cmap = ytcfg.get("yt", "default_colormap")
        grids = GridSource(data_source, alpha=alpha, cmap=cmap,
                            min_level=min_level, max_level=max_level)
        self.add_source(grids)
        return self

    def annotate_mesh_lines(self, color=None, alpha=1.0):
        """

        Modifies this Scene by drawing the mesh line boundaries
        on all MeshSources.

        Parameters
        ----------
        colors: array of ints, shape (4), optional
            The RGBA value to use to draw the mesh lines.
            Default is black.
        alpha : float, optional
            The opacity of the mesh lines. Default is 255 (solid).

        """
        for k, source in self.opaque_sources:
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
        colors: array-like, shape (3,4), optional
            The x, y, z RGBA values to use to draw the axes.
        alpha : float, optional
            The opacity of the vectors.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
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
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>>
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
        cast to arbitrary units using the ``input_units`` keyword argument.

        Parameters
        ----------

        input_array : iterable
            A tuple, list, or array to attach units to
        input_units : String unit specification, unit symbol object, or astropy
                      units object
            The units of the array. Powers must be specified using python syntax
            (cm**3, not cm^3).
        dtype : string or NumPy dtype object
            The dtype of the returned array data

        Examples
        --------

        >>> a = sc.arr([1, 2, 3], 'cm')
        >>> b = sc.arr([4, 5, 6], 'm')
        >>> a + b
        YTArray([ 401.,  502.,  603.]) cm
        >>> b + a
        YTArray([ 4.01,  5.02,  6.03]) m

        Arrays returned by this function know about the scene's unit system

        >>> a = sc.arr(np.ones(5), 'unitary')
        >>> a.in_units('Mpc')
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
        cast to arbitrary units using the ``input_units`` keyword argument.

        Parameters
        ----------

        input_scalar : an integer or floating point scalar
            The scalar to attach units to
        input_units : String unit specification, unit symbol object, or astropy
                      units
            The units of the quantity. Powers must be specified using python
            syntax (cm**3, not cm^3).
        dtype : string or NumPy dtype object
            The dtype of the array data.

        Examples
        --------

        >>> a = sc.quan(1, 'cm')
        >>> b = sc.quan(2, 'm')
        >>> a + b
        201.0 cm
        >>> b + a
        2.01 m

        Quantities created this way automatically know about the unit system
        of the scene

        >>> a = ds.quan(5, 'unitary')
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
        png = self._last_render.write_png(filename=None,
                                          sigma_clip=self._sigma_clip,
                                          background='black')
        self._sigma_clip = None
        return png

    def __repr__(self):
        disp = "<Scene Object>:"
        disp += "\nSources: \n"
        for k, v in iteritems(self.sources):
            disp += "    %s: %s\n" % (k, v)
        disp += "Camera: \n"
        disp += "    %s" % self.camera
        return disp

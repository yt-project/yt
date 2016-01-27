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


import numpy as np
from collections import OrderedDict
from yt.funcs import mylog, get_image_suffix
from yt.extern.six import iteritems, itervalues, string_types
from .camera import Camera
from .render_source import OpaqueSource, BoxSource, CoordinateVectorSource, \
    GridSource, RenderSource, MeshSource
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
    >>> cam = Camera(ds)
    >>> sc.camera = cam
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

    def __init__(self):
        r"""Create a new Scene instance"""
        super(Scene, self).__init__()
        self.sources = OrderedDict()
        self.camera = None
        # An image array containing the last rendered image of the scene
        self.last_render = None
        # A non-public attribute used to get around the fact that we can't
        # pass kwargs into _repr_png_()
        self._sigma_clip = None

    def get_source(self, source_num):
        """Returns the volume rendering source indexed by ``source_num``"""
        return list(itervalues(self.sources))[source_num]

    def _iter_opaque_sources(self):
        """
        Iterate over opaque RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in iteritems(self.sources):
            if isinstance(source, OpaqueSource) or \
                    issubclass(OpaqueSource, type(source)):
                yield k, source

    def _iter_transparent_sources(self):
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

        self.sources[keyname] = render_source

        return self

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
        self.last_render = bmp
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

        if self.last_render is None:
            self.render()

        mylog.info("Saving render %s", fname)
        self.last_render.write_png(fname, sigma_clip=sigma_clip)
 
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
        opaque = ZBuffer(empty, np.ones(empty.shape[:2]) * np.inf)

        for k, source in self._iter_opaque_sources():
            source.render(camera, zbuffer=opaque)
            im = source.zbuffer.rgba

        for k, source in self._iter_transparent_sources():
            im = source.render(camera, zbuffer=opaque)

        return im

    def camera():
        doc = r"""The camera property.

        This is the default camera that will be used when rendering. Can be set
        manually, but Camera type will be checked for validity.
        """

        def fget(self):
            cam = self._camera
            if cam is None:
                cam = Camera()
            self._camera = cam
            return self._camera

        def fset(self, value):
            # Should add better validation here
            self._camera = value

        def fdel(self):
            del self._camera
            self._camera = None
        return locals()
    camera = property(**camera())

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
                               color=None)
        self.add_source(box_source)
        return self

    def annotate_grids(self, data_source, alpha=0.3, cmap='algae',
                       min_level=None, max_level=None):
        r"""

        Modifies this scene by drawing the edges of the AMR grids.
        This adds a new BoxSource to the scene for each AMR grid 
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
        for k, source in self._iter_opaque_sources():
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
            self._sigma_clip = sigma_clip
            return self
        else:
            raise YTNotInsideNotebook

    def _repr_png_(self):
        if self.last_render is None:
            self.render()
        png = self.last_render.write_png(filename=None,
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

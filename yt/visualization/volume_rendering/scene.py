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


from yt.funcs import mylog
from camera import Camera
from render_source import OpaqueSource
from zbuffer_array import ZBuffer
import numpy as np


class SceneHandle(object):
    """docstring for SceneHandle"""
    def __init__(self, scene, camera, source, lens):
        mylog.debug("Entering %s" % str(self))
        self.scene = scene
        self.camera = camera
        self.source = source
        self.lens = lens

    def __repr__(self):
        desc = super(SceneHandle, self).__repr__()
        desc += str(self)
        return desc

    def __str__(self):
        desc = "Scene Handler\n"
        desc += ".scene: " + self.scene.__repr__() + "\n"
        desc += ".camera: " + self.camera.__repr__() + "\n"
        desc += ".source: " + self.source.__repr__() + "\n"
        desc += ".lens: " + self.lens.__repr__() + "\n"
        return desc


class Scene(object):

    """The Scene Class

    The Scene class is meant to be the primary container for the
    new volume rendering framework. A single scene may contain
    several Camera and RenderSource instances, and is the primary
    driver behind creating a volume rendering.

    """

    _current = None

    def __init__(self):
        r"""Create a new Scene instance.

        This sets up the basics needed to add sources and cameras.
        This does very little setup, and requires additional input
        to do anything useful.

        Parameters
        ----------
        None

        Examples
        --------
        >>> sc = Scene()

        """
        super(Scene, self).__init__()
        self.sources = {}
        self.camera = None

    def get_source(self, source_num):
        return self.sources.values()[source_num]

    def iter_opaque_sources(self):
        """
        Iterate over opaque RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.iteritems():
            if isinstance(source, OpaqueSource) or \
                    issubclass(OpaqueSource, type(source)):
                yield k, source

    def iter_transparent_sources(self):
        """
        Iterate over transparent RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.iteritems():
            if not isinstance(source, OpaqueSource):
                yield k, source

    def get_camera(self):
        """
        Use exisiting sources and their data sources to
        build a default camera. If no useful source is
        available, create a default Camera at 1,1,1 in the
        1,0,0 direction"""
        cam = self.camera
        if cam is None:
            cam = Camera()
        self.camera = cam
        return cam

    def add_source(self, render_source, keyname=None):
        """
        Add a render source to the scene.  This will autodetect the
        type of source.
        """
        if keyname is None:
            keyname = 'source_%02i' % len(self.sources)

        self.sources[keyname] = render_source

        return self

    def render(self, fname=None, clip_ratio=None, camera=None):
        r"""Render all sources in the Scene.

        Use the current state of the Scene object to render all sources
        currently in the scene.

        Parameters
        ----------
        fname: string, optional
            If specified, save the rendering as a bitmap to the file "fname".
            Default: None
        clip_ratio: float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()
        camera: :class:`Camera`, optional
            If specified, use a different :class:`Camera` to render the scene.

        Returns
        -------
        bmp: :class:`ImageArray`
            ImageArray instance of the current rendering image.

        Examples
        --------
        >>> sc = Scene()
        >>> # Add sources/camera/etc
        >>> im = sc.render('rendering.png')
        """
        if camera is None:
            camera = self.camera
        assert(camera is not None)
        self.validate()
        bmp = self.composite(camera=camera)
        if fname is not None:
            bmp.write_png(fname, clip_ratio=clip_ratio)
        return bmp

    def validate(self):
        for k, source in self.sources.iteritems():
            source.validate()
        return

    def composite(self, camera=None):
        if camera is None:
            camera = self.camera
        empty = camera.lens.new_image(camera)
        opaque = ZBuffer(empty, np.ones(empty.shape[:2]) * np.inf)

        for k, source in self.iter_opaque_sources():
            source.render(camera, zbuffer=opaque)

        for k, source in self.iter_transparent_sources():
            im = source.render(camera, zbuffer=opaque)
        return im

    def set_camera(self, camera):
        self.camera = camera

    def get_handle(self, key=None):
        """docstring for get_handle"""

        if key is None:
            key = self.sources.keys()[0]
        handle = SceneHandle(self, self.camera, self.sources[key],
                             self.sources[key].lens)
        return handle

    def __repr__(self):
        disp = "<Scene Object>:"
        disp += "\nSources: \n"
        for k, v in self.sources.iteritems():
            disp += "    %s: %s\n" % (k, v)
        disp += "Camera: \n"
        disp += "    %s" % self.camera
        return disp

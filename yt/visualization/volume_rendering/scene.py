"""


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


from yt.funcs import mylog
from yt.data_objects.static_output import Dataset
from camera import Camera
from render_source import VolumeSource, OpaqueSource
from yt.data_objects.api import ImageArray
from zbuffer_array import ZBuffer
import numpy as np


class SceneHandle(object):
    """docstring for SceneHandle"""
    def __init__(self, scene, camera, source, lens):
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

    """Skeleton Class for 3D Scenes"""

    _current = None

    def __init__(self):
        super(Scene, self).__init__()
        self.sources = {}
        self.camera = None

    def iter_opaque_sources(self):
        """
        Iterate over opaque RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.iteritems():
            if isinstance(source, OpaqueSource):
                yield k, source

    def iter_transparent_sources(self):
        """
        Iterate over opaque RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.iteritems():
            if not isinstance(source, OpaqueSource):
                yield k, source

    def get_default_camera(self):
        """
        Use exisiting sources and their data sources to
        build a default camera. If no useful source is
        available, create a default Camera at 1,1,1 in the
        1,0,0 direction"""
        for k, source in self.sources.iteritems():
            cam = source.get_default_camera()
            if cam is not None:
                break
        if cam is None:
            cam = Camera()
        self.default_camera = cam
        return cam

    def add_source(self, render_source, keyname=None):
        """
        Add a render source to the scene.  This will autodetect the
        type of source.
        """
        if keyname is None:
            keyname = 'source_%02i' % len(self.sources)

        render_source.set_scene(self)

        self.sources[keyname] = render_source

        return self

    def render(self, fname=None, clip_ratio=None):
        self.validate()
        ims = {}
        for k, v in self.sources.iteritems():
            v.validate()
            print 'Running', k, v
            ims[k] = v.render()

        bmp = np.zeros_like(ims.values()[0])
        for k, v in ims.iteritems():
            bmp += v
        bmp = ImageArray(bmp.d)
        assert(isinstance(bmp, ImageArray))

        if fname is not None:
            bmp.write_png(fname, clip_ratio=clip_ratio)
        return bmp

    def validate(self):
        for k, source in self.sources.iteritems():
            source.validate()
        return

    def composite(self):
        opaque = ZBuffer(
            np.zeros(self.camera.resolution[0],
                     self.camera.resolution[1],
                     4),
            np.ones(self.camera.resolution) * np.inf)

        for k, source in self.iter_opaque_sources():
            if source.zbuffer is not None:
                opaque = opaque + source.zbuffer

        for k, source in self.iter_transparent_sources():
            source.render(zbuffer=opaque)
            opaque = opaque + source.zbuffer
        pass

    def set_default_camera(self, camera):
        self.default_camera = camera

    def get_handle(self, key=None):
        """docstring for get_handle"""

        if key is None:
            key = self.sources.keys()[0]
        handle = SceneHandle(self, self.camera, self.sources[key],
                             self.sources[key].lens)
        return handle


def create_volume_rendering(data_source, field=None):
    if isinstance(data_source, Dataset):
        pf = data_source
        data_source = data_source.all_data()
    else:
        pf = data_source.pf

    sc = Scene()
    camera = Camera(data_source)
    if field is None:
        pf.field_list
        field = pf.field_list[0]
        mylog.info('Setting default field to %s' % field.__repr__())
    render_source = VolumeSource(data_source, field)

    sc.set_camera(camera)
    sc.add_source(render_source)
    render_source.build_defaults()
    return sc

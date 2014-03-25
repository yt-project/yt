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
import numpy as np

class SceneHandle(object):
    """docstring for SceneHandle"""
    def __init__(self, scene, camera, source, engine):
        self.scene = scene
        self.camera = camera
        self.source = source
        self.engine = engine

    def __repr__(self):
        desc = super(SceneHandle, self).__repr__()
        desc += str(self)
        return desc

    def __str__(self):
        desc = "Scene Handler\n"
        desc += ".scene: " + self.scene.__repr__() + "\n"
        desc += ".camera: " + self.camera.__repr__() + "\n"
        desc += ".source: " + self.source.__repr__() + "\n"
        desc += ".engine: " + self.engine.__repr__() + "\n"
        return desc

class Scene(object):

    """Skeleton Class for 3D Scenes"""

    _current = None

    def __init__(self):
        super(Scene, self).__init__()
        self.camera = None
        self.sources = {}
        self.camera_path = None

    def set_camera(self, camera):
        self.camera = camera

        for source in self.sources.values():
            source.set_camera(self.camera)

    def get_handle(self, key=None):
        """docstring for get_handle"""

        if key is None:
            key = self.sources.keys()[0]
        handle = SceneHandle(self, self.camera, self.sources[key],
                             self.sources[key].engine)
        return handle

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


    def validate(self):
        if self.camera is None:
            for k, source in self.sources.iteritems():
                try:
                    self.camera = Camera(source.data_source)
                    return
                except:
                    pass
                raise RuntimeError("Couldn't build default camera")
        return

    def request(self):
        pass

    def composite(self):
        pass

    @property
    def current(self):
        if self._current is None:
            self.request()
        return self._current

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

class RenderScene(Scene):

    """docstring for RenderScene"""

    def __init__(self, data_source, field=None):
        super(RenderScene, self).__init__()
        if isinstance(data_source, Dataset):
            self.ds = data_source
            data_source = data_source.all_data()
        else:
            self.ds = data_source.pf

        self.data_source = data_source
        self.camera = Camera(data_source)
        self.field = field
        self.render_sources = {}
        self.default_setup()

    def default_setup(self):
        """docstring for default_setup"""
        if self.field is None:
            self.ds.field_list
            self.field = self.ds.field_list[0]
            mylog.info('Setting default field to %s' % self.field.__repr__())

        if self.data_source:
            render_source = VolumeSource(self.data_source, self.field)
            self.add_source(render_source)
            render_source.build_defaults()

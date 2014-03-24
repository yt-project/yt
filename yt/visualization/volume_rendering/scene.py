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


class Scene(object):

    """Skeleton Class for 3D Scenes"""

    _current = None

    def __init__(self):
        super(Scene, self).__init__()
        self.datasets = []
        self.camera = None
        self.sources = {}
        self.camera_path = None

    def set_camera(self, camera):
        self.camera = camera

        for source in self.sources.values():
            source.set_camera(self.camera)

    def setup_camera_links(self):
        """
        The camera object needs to be linked to:
            * Engines
            * Render Sources
        """
        if self.camera is None:
            raise RuntimeError("Camera does not exist")

        for source in self.sources.values():
            source.set_camera(self.camera)

    def iter_opaque_sources(self):
        """
        Iterate over opaque RenderSource objects,
        returning a tuple of (key, source)
        """
        for k, source in self.sources.iteritems():
            if isinstance(source, OpaqueSource):
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

    def register_dataset(self, ds):
        """Add a dataset to the scene"""
        self.datasets.append(ds)

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

    def render(self):
        self.validate()
        ims = {}
        for k, v in self.sources.iteritems():
            v.validate()
            print 'Running', k, v
            ims[k] = v.request()

        return ims


class RenderScene(Scene):

    """docstring for RenderScene"""

    def __init__(self, data_source=None, field=None):
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



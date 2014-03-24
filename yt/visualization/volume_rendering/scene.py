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
from render_source import VolumeSource


class Scene(object):

    """Skeleton Class for 3D Scenes"""

    _current = None

    def __init__(self):
        super(Scene, self).__init__()
        self.datasets = []
        self.camera = None
        self.sources = {}
        self.camera_path = None

    def request(self):
        pass

    @property
    def current(self):
        if self._current is None:
            self.request()
        return self._current

    def register_dataset(self, ds):
        """Add a dataset to the scene"""
        self.datasets.append(ds)

    def add_volume_rendering(self):
        """docstring for add_volume_rendering"""
        pass

    def add_slice(self):
        """docstring for add_slice"""
        pass

    def add_streamlines(self):
        """docstring for add_streamlines"""
        pass


class RenderScene(Scene):

    """docstring for RenderScene"""

    def __init__(self, data_source=None, field=None):
        super(RenderScene, self).__init__()
        if isinstance(data_source, Dataset):
            self.ds = data_source
            data_source = data_source.all_data()
        else:
            self.ds = data_source.pf

        print 'DATA SOURCE: ', data_source
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
            print 'WHAT FIELD AM I: ', self.field
            mylog.info('Setting default field to %s' % self.field.__repr__())

        if self.data_source:
            self.render_sources['vr1'] = \
                VolumeSource(self, self.data_source, self.field)

    def render(self):
        ims = {}
        for k, v in self.render_sources.iteritems():
            print 'Running', k, v
            ims[k] = v.request()

        return ims

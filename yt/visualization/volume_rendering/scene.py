"""


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


from yt.data_objects.static_output import Dataset


class Scene(object):

    """Skeleton Class for 3D Scenes"""

    _current = None

    def __init__(self):
        super(Scene, self).__init__()
        self.datasets = []
        self.camera = None
        self.sources = {}

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

    def __init__(self, data_source=None,):
        super(RenderScene, self).__init__()
        if isinstance(data_source, Dataset):
            self.ds = data_source
            data_source = ds.all_data()
        else:
            self.ds = data_source.pf
        self.data_source = data_source

        self.camera = Camera(data_source)
        self.render_sources = {}
        self.camera_path = CameraPath()

        if self.data_source:
            self.render_sources['vr1'] = VolumeSource(self.data_source)




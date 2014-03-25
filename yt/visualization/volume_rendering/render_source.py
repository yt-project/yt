"""
RenderSource Class



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.data_objects.api import ImageArray
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree
from transfer_function_helper import TransferFunctionHelper
from engine import PlaneParallelEngine
from camera import Camera


class RenderSource(ParallelAnalysisInterface):

    """Base Class for Render Sources. Will be inherited for volumes,
       streamlines, etc"""

    def __init__(self):
        super(RenderSource, self).__init__()
        self.opaque = False

    def request(self, *args, **kwargs):
        """returns a new ImageArray"""
        pass

    def setup(self):
        """Set up data needed to render"""
        pass


class OpaqueSource(RenderSource):
    """docstring for OpaqueSource"""
    def __init__(self):
        super(OpaqueSource, self).__init__()
        self.opaque = True


class VolumeSource(RenderSource):

    """docstring for VolumeSource"""

    def __init__(self, data_source, field=None):
        super(VolumeSource, self).__init__()
        self.data_source = data_source
        self.field = field
        self.scene = None
        self.volume = None
        self.current_image = None
        self.engine = None

        # In the future these will merge
        self.transfer_function = None
        self.tfh = None
        self.build_default_volume()

    def build_defaults(self):
        if self.data_source is not None:
            self.build_default_transfer_function()
            self.build_default_engine()

    def validate(self):
        """Make sure that all dependencies have been met"""
        if self.scene is None:
            raise RuntimeError("Scene not initialized")

        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

        if self.volume is None:
            raise RuntimeError("Volume not initialized")

        if self.engine is None:
            raise RuntimeError("Engine not initialized")

        if self.transfer_function is None:
            raise RuntimeError("Transfer Function not Supplied")
        self.setup()

    def build_default_transfer_function(self):
        self.tfh = \
            TransferFunctionHelper(self.data_source.pf)
        self.tfh.set_field(self.field)
        self.tfh.build_transfer_function()
        self.tfh.setup_default()
        self.transfer_function = self.tfh.tf

    def build_default_engine(self):
        self.engine = PlaneParallelEngine(self.scene, self)

    def build_default_volume(self):
        self.volume = AMRKDTree(self.data_source.pf,
                                data_source=self.data_source)
        log_fields = [self.data_source.pf.field_info[self.field].take_log]
        self.volume.set_fields([self.field], log_fields, True)

    def set_scene(self, scene):
        self.scene = scene
        if self.engine is not None:
            self.engine.set_camera(scene.camera)

    def set_camera(self, camera):
        """Set camera in this object, as well as any attributes"""
        self.engine.set_camera(camera)

    def prepare(self):
        """prepare for rendering"""
        self.scene.validate()
        self.new_image()

    def new_image(self):
        cam = self.scene.camera
        if cam is None:
            cam = Camera(self.data_source)
            self.scene.camera = cam
        self.current_image = ImageArray(
            np.zeros((cam.resolution[0], cam.resolution[1],
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def teardown(self):
        """docstring for teardown"""
        pass

    def add_sampler(self, sampler):
        """docstring for add_sampler"""
        pass

    def request(self, *args, **kwargs):
        """docstring for request"""
        self.prepare()
        self.engine.run()
        return self.current_image

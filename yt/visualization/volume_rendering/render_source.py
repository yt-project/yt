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
from yt.funcs import mylog
from yt.data_objects.api import ImageArray
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree
from transfer_function_helper import TransferFunctionHelper
from lens import PlaneParallelLens
from camera import Camera


class RenderSource(ParallelAnalysisInterface):

    """Base Class for Render Sources. Will be inherited for volumes,
       streamlines, etc"""

    def __init__(self):
        super(RenderSource, self).__init__()
        self.opaque = False
        self.lens = None
        self.zbuffer = None

    def setup(self):
        """Set up data needed to render"""
        pass

    def set_scene(self, scene):
        self.scene = scene
        if self.lens is not None:
            self.lens.set_camera(scene.camera)

    def render(self, zbuffer=None):
        pass

    def validate(self):
        pass

    def new_image(self):
        pass

    def prepare(self):
        pass

    def get_default_camera(self):
        """If possible, create a camera based on the render source"""
        return None


class OpaqueSource(RenderSource):
    """docstring for OpaqueSource"""
    def __init__(self):
        super(OpaqueSource, self).__init__()
        self.opaque = True

    def set_zbuffer(self, zbuffer):
        self.zbuffer = zbuffer

class VolumeSource(RenderSource):

    """docstring for VolumeSource"""

    def __init__(self, data_source, field=None):
        super(VolumeSource, self).__init__()
        self.data_source = data_source
        self.field = field
        self.scene = None
        self.volume = None
        self.current_image = None
        self.lens = None

        # In the future these will merge
        self.transfer_function = None
        self.tfh = None
        self.build_default_volume()

    def build_defaults(self):
        if self.data_source is not None:
            self.build_default_transfer_function()

    def validate(self):
        """Make sure that all dependencies have been met"""
        if self.scene is None:
            raise RuntimeError("Scene not initialized")

        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

        if self.volume is None:
            raise RuntimeError("Volume not initialized")

        #if self.lens is None:
        #    raise RuntimeError("Lens not initialized")

        if self.transfer_function is None:
            raise RuntimeError("Transfer Function not Supplied")
        self.setup()

    def switch_field(self, field):
        self.volume.set_fields([self.field], log_fields, True)


    def build_default_transfer_function(self):
        self.tfh = \
            TransferFunctionHelper(self.data_source.pf)
        self.tfh.set_field(self.field)
        self.tfh.build_transfer_function()
        self.tfh.setup_default()
        self.transfer_function = self.tfh.tf

    def prepare(self):
        """prepare for rendering"""
        self.scene.validate()
        self.new_image()

    def build_default_volume(self):
        self.volume = AMRKDTree(self.data_source.pf,
                                data_source=self.data_source)
        log_fields = [self.data_source.pf.field_info[self.field].take_log]
        mylog.debug('Log Fields:' + str(log_fields))
        self.volume.set_fields([self.field], log_fields, True)

    def set_camera(self, camera):
        """Set camera in this object, as well as any attributes"""
        self.lens.set_camera(camera)

    def teardown(self):
        """docstring for teardown"""
        pass

    def add_sampler(self, sampler):
        """docstring for add_sampler"""
        pass

    def render(self, zbuffer=None):
        """docstring for request"""
        self.zbuffer = zbuffer
        self.prepare()
        self.lens.run()

        self.camera_updated()
        total_cells = 0
        if self.double_check:
            for brick in self.render_source.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        view_pos = self.front_center + \
            self.camera.unit_vectors[2] * 1.0e6 * self.camera.width[2]
        for brick in self.render_source.volume.traverse(view_pos):
            self.sampler(brick, num_threads=self.num_threads)
            total_cells += np.prod(brick.my_data[0].shape)

        self.current_image = \
            self.finalize_image(self.sampler.aimage)

        return self.current_image


"""
Engine Classes



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.lib.grid_traversal import \
    VolumeRenderSampler
from yt.units.yt_array import YTArray
import numpy as np


class Engine(ParallelAnalysisInterface):

    """docstring for Engine"""

    def __init__(self, ):
        super(Engine, self).__init__()

    def camera_updated(self):
        """docstring for update_camera"""
        pass


class PlaneParallelEngine(Engine):

    """docstring for PlaneParallelEngine"""

    def __init__(self, scene, render_source):
        super(PlaneParallelEngine, self).__init__()
        self.scene = scene
        self.camera = scene.camera
        self.render_source = render_source
        self.transfer_function = self.render_source.transfer_function
        self.sub_samples = 5
        self.num_threads = 1
        self.double_check = False
        self.box_vectors = None
        self.origin = None
        self.back_center = None
        self.front_center = None

        if scene.camera:
            self._setup_box_properties()
        self.sampler = self.get_sampler()

    def set_camera(self, camera):
        """set the camera for this engine"""
        self.camera = camera

    def camera_updated(self):
        if self.camera._moved:
            self._setup_box_properties()
            self.sampler = self.get_sampler()
            self.camera._moved = False

    def _setup_box_properties(self):
        unit_vectors = self.camera.unit_vectors
        width = self.camera.width
        center = self.camera.focus
        self.box_vectors = YTArray([unit_vectors[0] * width[0],
                                    unit_vectors[1] * width[1],
                                    unit_vectors[2] * width[2]])
        self.origin = center - 0.5 * width.dot(YTArray(unit_vectors, ""))
        self.back_center = center - 0.5 * width[2] * unit_vectors[2]
        self.front_center = center + 0.5 * width[2] * unit_vectors[2]
        mylog.debug('Setting box properties')
        mylog.debug(self.back_center)
        mylog.debug(self.front_center)

    def get_sampler(self):
        self.render_source.prepare()
        image = self.render_source.current_image
        rotp = np.concatenate([self.scene.camera.inv_mat.ravel('F'),
                               self.back_center.ravel()])
        args = (rotp, self.box_vectors[2], self.back_center,
                (-self.camera.width[0] / 2.0, self.camera.width[0] / 2.0,
                 -self.camera.width[1] / 2.0, self.camera.width[1] / 2.0),
                image, self.camera.unit_vectors[
                    0], self.camera.unit_vectors[1],
                np.array(self.camera.width, dtype='float64'),
                self.render_source.transfer_function, self.sub_samples)
        sampler = VolumeRenderSampler(*args)
        return sampler

    def run(self):
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

        self.render_source.current_image = \
            self.finalize_image(self.sampler.aimage)
        return

    def finalize_image(self, image):
        cam = self.scene.camera
        view_pos = self.front_center + cam.unit_vectors[2] * \
            1.0e6 * cam.width[2]
        image = self.render_source.volume.reduce_tree_images(image, view_pos)
        if self.transfer_function.grey_opacity is False:
            image[:, :, 3] = 1.0
        return image

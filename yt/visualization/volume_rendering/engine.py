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

from yt.funcs import get_pbar
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.lib.grid_traversal import \
    VolumeRenderSampler
import numpy as np


class Engine(ParallelAnalysisInterface):

    """docstring for Engine"""

    def __init__(self, ):
        super(Engine, self).__init__()


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
        self.sampler = self.get_sampler()

    def get_sampler(self):
        image = self.render_source.current_image
        rotp = np.concatenate([self.scene.camera.inv_mat.ravel('F'),
                               self.scene.camera.back_center.ravel()])
        args = (rotp, self.camera.box_vectors[2], self.camera.back_center,
                (-self.camera.width[0] / 2.0, self.camera.width[0] / 2.0,
                 -self.camera.width[1] / 2.0, self.camera.width[1] / 2.0),
                image, self.camera.unit_vectors[
                    0], self.camera.unit_vectors[1],
                np.array(self.camera.width, dtype='float64'),
                self.transfer_function, self.sub_samples)
        sampler = VolumeRenderSampler(*args)
        return sampler

    def run(self):
        total_cells = 0
        if self.double_check:
            for brick in self.render_source.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        view_pos = self.camera.front_center + \
            self.camera.unit_vectors[2] * 1.0e6 * self.camera.width[2]
        for brick in self.render_source.volume.traverse(view_pos):
            self.sampler(brick, num_threads=self.num_threads)
            total_cells += np.prod(brick.my_data[0].shape)

        self.render_source.current_image = \
            self.finalize_image(self.sampler.aimage)
        return

    def finalize_image(self, image):
        cam = self.scene.camera
        view_pos = cam.front_center + cam.unit_vectors[2] * \
            1.0e6 * cam.width[2]
        image = self.render_source.volume.reduce_tree_images(image, view_pos)
        if self.transfer_function.grey_opacity is False:
            image[:, :, 3] = 1.0
        return image

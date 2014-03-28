"""
Lens Classes



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
from camera import Camera
from yt.units.yt_array import YTArray
import numpy as np


class Lens(ParallelAnalysisInterface):

    """docstring for Lens"""

    def __init__(self, ):
        super(Lens, self).__init__()

    def camera_updated(self):
        """docstring for update_camera"""
        pass


class PlaneParallelLens(Lens):

    """docstring for PlaneParallelLens"""

    def __init__(self, scene, render_source):
        super(PlaneParallelLens, self).__init__()
        self.scene = scene
        self.camera = scene.camera
        self.render_source = render_source
        self.sub_samples = 5
        self.num_threads = 0
        self.double_check = False
        self.box_vectors = None
        self.origin = None
        self.back_center = None
        self.front_center = None

        if scene.camera:
            self._setup_box_properties()
        self.sampler = self.get_sampler()

    def set_camera(self, camera):
        """set the camera for this lens"""
        self.camera = camera

    def camera_updated(self):
        if self.camera._moved:
            self._setup_box_properties()
            self.sampler = self.get_sampler()
            self.camera._moved = False

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
        self._setup_box_properties()
        kwargs = {}
        if self.render_source.zbuffer is not None:
            kwargs['zbuffer'] = self.render_source.zbuffer.z
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
        sampler = VolumeRenderSampler(*args, **kwargs)
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



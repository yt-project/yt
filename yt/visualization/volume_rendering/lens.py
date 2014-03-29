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
from yt.data_objects.image_array import ImageArray
import numpy as np


class Lens(ParallelAnalysisInterface):

    """docstring for Lens"""

    def __init__(self, ):
        super(Lens, self).__init__()


class PlaneParallelLens(Lens):

    """docstring for PlaneParallelLens"""

    def __init__(self):
        super(PlaneParallelLens, self).__init__()
        self.sub_samples = 5
        self.num_threads = 0
        self.double_check = False
        self.box_vectors = None
        self.origin = None
        self.back_center = None
        self.front_center = None
        self.sampler = None

    def expose(self, scene, camera, render_source):
        self.setup_box_properties(camera)
        self.sampler = self.get_sampler(scene, camera, render_source)
        self.cast_rays(camera, self.sampler, render_source)

    def new_image(self, camera):
        cam = camera
        if cam is None:
            cam = Camera(self.data_source)
        self.current_image = ImageArray(
            np.zeros((cam.resolution[0], cam.resolution[1],
                      4), dtype='float64', order='C'),
            info={'imtype': 'rendering'})
        return self.current_image

    def setup_box_properties(self, camera):
        unit_vectors = camera.unit_vectors
        width = camera.width
        center = camera.focus
        self.box_vectors = YTArray([unit_vectors[0] * width[0],
                                    unit_vectors[1] * width[1],
                                    unit_vectors[2] * width[2]])
        self.origin = center - 0.5 * width.dot(YTArray(unit_vectors, ""))
        self.back_center = center - 0.5 * width[2] * unit_vectors[2]
        self.front_center = center + 0.5 * width[2] * unit_vectors[2]

    def get_sampler(self, scene, camera, render_source):
        kwargs = {}
        if render_source.zbuffer is not None:
            kwargs['zbuffer'] = render_source.zbuffer.z
        render_source.prepare()
        image = render_source.current_image
        image = self.new_image(camera)
        rotp = np.concatenate([camera.inv_mat.ravel('F'),
                               self.back_center.ravel()])
        args = (rotp, self.box_vectors[2], self.back_center,
                (-camera.width[0] / 2.0, camera.width[0] / 2.0,
                 -camera.width[1] / 2.0, camera.width[1] / 2.0),
                image, camera.unit_vectors[
                    0], camera.unit_vectors[1],
                np.array(camera.width, dtype='float64'),
                render_source.transfer_function, self.sub_samples)
        sampler = VolumeRenderSampler(*args, **kwargs)
        return sampler

    def cast_rays(self, camera, sampler, render_source):
        mylog.debug("Casting rays")
        total_cells = 0
        if self.double_check:
            for brick in render_source.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        # This is a hack that should be replaced by an alternate plane-parallel
        # traversal. Put the camera really far away so that the effective
        # viewpoint is infinitely far away, making for parallel rays.
        view_pos = self.front_center + \
            camera.unit_vectors[2] * 1.0e6 * camera.width[2]

        for brick in render_source.volume.traverse(view_pos):
            sampler(brick, num_threads=self.num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
        mylog.debug("Done casting rays")

        render_source.current_image = \
            self.finalize_image(camera, render_source,
                                         self.sampler.aimage)
        return

    def finalize_image(self, camera, render_source, image):
        cam = camera
        view_pos = self.front_center + cam.unit_vectors[2] * \
            1.0e6 * cam.width[2]
        image = render_source.volume.reduce_tree_images(image, view_pos)
        if render_source.transfer_function.grey_opacity is False:
            image[:, :, 3] = 1.0
        return image

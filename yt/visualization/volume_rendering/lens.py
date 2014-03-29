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
from yt.units.yt_array import YTArray
from yt.data_objects.image_array import ImageArray
import numpy as np


class Lens(ParallelAnalysisInterface):

    """docstring for Lens"""

    def __init__(self, ):
        mylog.debug("Entering %s" % str(self))
        super(Lens, self).__init__()
        self.viewpoint = None


class PlaneParallelLens(Lens):

    """docstring for PlaneParallelLens"""

    def __init__(self):
        mylog.debug("Entering %s" % str(self))
        super(PlaneParallelLens, self).__init__()
        self.sub_samples = 5
        self.num_threads = 0
        self.double_check = False
        self.box_vectors = None
        self.origin = None
        self.back_center = None
        self.front_center = None
        self.sampler = None
        self.viewpoint = None

    def set_camera(self, camera):
        self.setup_box_properties(camera)

    def expose(self, scene, camera, render_source):
        self.setup_box_properties(camera)
        self.sampler = self.get_sampler(scene, camera, render_source)
        self.cast_rays(camera, self.sampler, render_source)

    def new_image(self, camera):
        self.current_image = ImageArray(
            np.zeros((camera.resolution[0], camera.resolution[1],
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

        # This is a hack that should be replaced by an alternate plane-parallel
        # traversal. Put the camera really far away so that the effective
        # viewpoint is infinitely far away, making for parallel rays.
        self.viewpoint = self.front_center + \
            camera.unit_vectors[2] * 1.0e6 * camera.width[2]

    def get_sampler_params(self, camera):
        sampler_params =\
            dict(vp_pos=np.concatenate([camera.inv_mat.ravel('F'),
                                        self.back_center.ravel()]),
                 vp_dir=self.box_vectors[2],  # All the same
                 center=self.back_center,
                 bounds=(-camera.width[0] / 2.0, camera.width[0] / 2.0,
                         -camera.width[1] / 2.0, camera.width[1] / 2.0),
                 x_vec=camera.unit_vectors[0],
                 y_vec=camera.unit_vectors[1],
                 width=np.array(camera.width, dtype='float64'),
                 image=self.new_image(camera))
        return sampler_params


class PerspectiveLens(Lens):
    """docstring for PerspectiveLens"""
    def __init__(self):
        super(PerspectiveLens, self).__init__()
        raise NotImplementedError


class FisheyeLens(Lens):
    """docstring for FisheyeLens"""
    def __init__(self):
        super(FisheyeLens, self).__init__()
        raise NotImplementedError

lenses = {'plane-parallel': PlaneParallelLens,
          'perspective': PerspectiveLens,
          'fisheye': FisheyeLens}

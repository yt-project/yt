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

from yt.data_objects.api import ImageArray
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree


class RenderSource(ParallelAnalysisInterface):

    """Base Class for Render Sources. Will be inherited for volumes,
       streamlines, etc"""

    def __init__(self):
        super(RenderSource, self).__init__()
        pass

    def request(self, *args, **kwargs):
        """returns a new ImageArray"""
        pass


class VolumeSource(RenderSource):

    """docstring for VolumeSource"""

    def __init__(self, data_source):
        super(VolumeSource, self).__init__()
        self.data_source = data_source
        self.volume = None
        self.setup()

    def validate(self):
        """docstring for validate"""
        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

        if self.volume is None:
            raise RuntimeError("Volume not initialized")

    def setup(self):
        """setup VolumeSource"""
        self.volume = AMRKDTree(self.data_source)

    def teardown(self):
        """docstring for teardown"""
        pass

    def add_sampler(self, sampler):
        """docstring for add_sampler"""
        pass

    def request(self, *args, **kwargs):
        """docstring for request"""
        pass

    def _render(self, double_check, num_threads, image, sampler):
        pbar = get_pbar("Ray casting",
                        (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        for brick in self.volume.traverse(self.front_center):
            sampler(brick, num_threads=num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
            pbar.update(total_cells)

        pbar.finish()
        image = self.finalize_image(sampler.aimage)
        return image



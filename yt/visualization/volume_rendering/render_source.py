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
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree
from transfer_function_helper import TransferFunctionHelper
from engine import PlaneParallelEngine


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

    def __init__(self, scene, data_source, field=None):
        super(VolumeSource, self).__init__()
        self.scene = scene
        self.data_source = data_source
        self.volume = None
        self.field = field
        self.current_image = None
        self.engine = None
        self.setup()

        # In the future these will merge
        self.transfer_function = None
        self.tfh = None
        self.setup()

    def validate(self):
        """docstring for validate"""
        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

        if self.volume is None:
            raise RuntimeError("Volume not initialized")

        if self.transfer_function is None:
            raise RuntimeError("Transfer Function not Supplied")

    def setup(self):
        """setup VolumeSource"""
        self.current_image = self.new_image()
        self.tfh = \
            TransferFunctionHelper(self.data_source.pf)
        self.tfh.set_field(self.field)
        self.tfh.set_bounds([0.0, 1.0])
        self.tfh.set_log(False)
        self.tfh.build_transfer_function()
        self.tfh.setup_default()
        self.transfer_function = self.tfh.tf
        self.engine = PlaneParallelEngine(self.scene, self)
        self.volume = AMRKDTree(self.data_source.pf,
                                data_source=self.data_source)
        self.volume.initialize_source([self.field], [False], True)

    def new_image(self):
        cam = self.scene.camera
        image = np.zeros((cam.resolution[0],
                          cam.resolution[1], 4),
                         dtype='float64', order='C')
        return image

    def teardown(self):
        """docstring for teardown"""
        pass

    def add_sampler(self, sampler):
        """docstring for add_sampler"""
        pass

    def request(self, *args, **kwargs):
        """docstring for request"""
        self.engine.run()
        return self.current_image

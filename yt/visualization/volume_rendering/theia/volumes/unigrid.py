#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import pycuda.driver as drv
import yt.visualization.volume_rendering.theia.helpers.cuda as cdh

from yt.visualization.volume_rendering.theia.volumes.volume import Volume

import numpy as np

class Unigrid(Volume):
	def __init__(self, array = None, allocate = False,
			     flags = drv.TRSF_NORMALIZED_COORDINATES,
			     filter_mode = drv.filter_mode.LINEAR,
			     address_mode = drv.address_mode.CLAMP):
            Volume.__init__(self)

            self.flags = flags
            self.filter_mode = filter_mode
            self.address_mode = address_mode

            self.volume = None
            self.cuda_volume_array = None

            if array != None:
                  self.set_array(array)
                  if allocate == True:
                        self.allocate_and_copy_to_gpu()


	def set_filter_mode(self, mode):
		self.filter_mode = mode

	def set_address_mode(self, mode):
		self.address_mode = mode

	def set_flags(self, flags):
		self.flags = flags

	def set_array(self, np_array):
		self.volume = np_array

	def allocate_and_copy_to_gpu(self):
            if self.cuda_volume_array:
                  self.deallocate_from_gpu()

            self.cuda_volume_array = cdh.np3d_to_device_array(self.volume)

        def __del__(self) :
            self.deallocate_from_gpu()

	def deallocate_from_gpu(self):
		self.cuda_volume_array.free()
		self.cuda_volume_array = None

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


"""Linear Transfer Function
This module defines a linear-type transfer function.  This type
of transfer function is constructed from a contiguous array of
RGBA values.  In this case, the array input should represent a
numpy array.
"""

import pycuda.driver as drv
import yt.visualization.volume_rendering.theia.transfer.helper as tfh

import numpy as np

class LinearTransferFunction:
	def __init__(self, arrays = None,
	             filter_mode = drv.filter_mode.LINEAR,
				 address_mode = drv.address_mode.CLAMP,
				 flags = drv.TRSF_NORMALIZED_COORDINATES,
				 range = (0.0, 1.0)):

                (low, up) = range
                self.scale  = 1.0/(up - low)
                self.offset = low
                if (low < 0 and up < 0):
                    self.scale  *= -1.0


		self.range = range
		self.flags = flags
		self.filter_mode = filter_mode
		self.address_mode = address_mode
		self.interleaved = None

		if arrays != None:
			self.arrays_to_transfer_function(arrays)


	def set_filter_mode(self, mode):
		"""Sets the Texture Filter Mode

		"""
		self.filter_mode = mode

	def set_address_mode(self, mode):
		self.address_mode = mode

	def set_flags(self, flags):
		self.flags = flags

	def set_range(self, range):
		self.range = range

		self.offset = tfh.bounds_to_offset(range)
		self.scale = tfh.bounds_to_scale(range)


	def arrays_to_transfer_function(self, arrays):

		(r_array, g_array, b_array, a_array) = arrays
		(size, ) = r_array.shape


		self.interleaved = np.asarray(np.dstack((r_array.reshape(1, size),
										g_array.reshape(1, size),
										b_array.reshape(1, size),
										a_array.reshape(1, size))),
                                                            dtype=np.float32, order='F')
		self.cuda_transfer_array = drv.make_multichannel_2d_array(np.asarray(self.interleaved.transpose((2,1,0)),
													   dtype=np.float32, order='F'), order='F')

	def yt_to_function(self, tf):
		size = tf.nbins
		self.interleaved = np.asarray(np.dstack((tf.red.y.reshape(1, size),
					         	 tf.green.y.reshape(1, size),
							 tf.blue.y.reshape(1, size),
							 tf.alpha.y.reshape(1, size))),
					      dtype=np.float32, order='F')

		self.cuda_transfer_array = drv.make_multichannel_2d_array(np.asarray(self.interleaved.transpose((2,1,0)),
			 							     dtype=np.float32, order='F'), order='F')

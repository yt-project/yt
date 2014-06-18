#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


class Surface:
        """
         Containe for holding image array, results of volume rendering.

         Parameters
         ----------
         device_ptr: uint64
                     A pycuda allocation object representing a point to device memory.
         bounds:     
                     A tuple representing the shape of the underlying data, often
                       numpy.shape.
        """
	def __init__(self):
		self.device_ptr = None
		self.bounds = None

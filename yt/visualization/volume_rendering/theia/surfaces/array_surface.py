#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import pycuda.driver as drv

import numpy as np

from   yt.visualization.volume_rendering.theia.surfaces.surface import Surface


"""
An array surface is a type of surface object that is bound to a specified
n-dimensional numpy array.	It contains functions that allow data to be
readily copied back and forth between host and device.
"""

class ArraySurface(Surface):
	def __init__(self, array=None, size=None, copy=True) :
            """
                Creates a new ArraySurface.	 May optionally be given initial values.

                Parameters
                ----------
                array:      An ndarray to initialize local_array to. Will copy to
                          device by default.

                size:	  A shape that would serve as a basis for an ndarray.
                          this argument is not necessary if array argument is specified
                          as the shape can be inferred from it.

                copy:	  When True, will copy local array to device.
            """
            Surface.__init__(self)

            self.local_array = None
            self.bounds      = size

            if (array != None):
                self.set_array(array, copy=copy)
            else:
                if size != None:
                      self.set_array(np.zeros(size, dtype=np.uint32), copy=copy)


	def device_to_local(self):
		"""
		    Copies device memory to local host memory.
		"""

		self.local_array = drv.from_device_like(self.device_ptr, self.local_array)


	def local_to_device(self, early_free=True):
		"""
		    Copies locally set array data to device memory.
                    Parameters
                    ----------
		    early_free: boolean 
                            When True, will automatically free all device memory
		            before allocation, preventing double-allocation for a short time.
		"""

		if self.device_ptr != None and early_free == True:
			self.device_ptr.free()

		self.device_ptr = drv.to_device(self.local_array)


	def set_array(self, array, copy=True):
            """
                Sets the local array in the Surface to the first argument, array, a
                a numpy array.
                Parameters
                ----------
                array : 2D numpy array
                     Container defining the surface.
                copy: boolean
                     When True, will immediately copy local data to device.
            """
            self.local_array = drv.pagelocked_empty_like(array)
            self.bounds = array.shape
            if copy == True:
                  self.local_to_device()


	def get_array(self, copy=True):
		"""
		     Gets a copy of the local array and returns it.

                     Parameters
                     ----------
 
		     copy: boolean
                           When True, will immediately copy device data
                           into local memory.
		"""

		if copy == True:
			self.device_to_local()
		return self.local_array


	def get_bounds(self):
		"""
		    Returns the shape of the local ndarray.
		"""
		return self.bounds

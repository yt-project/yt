#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import pycuda.driver as drv
import pyRGBA.theia.helpers.cuda as cdh
import numpy as np

from pyRGBA.theia.transfer.linear_transfer import LinearTransferFunction

from pyRGBA.theia.volumes.cube import Unigrid
from pyRGBA.theia.surfaces.array_surface import ArraySurface
import os

class FrontToBackRaycaster:
    r"""
    This is a basic camera intended to be a base class. 
    The camera provies utilities for controling the view point onto the data.


    Parameters
    ----------
    pos   : numpy array

    Example:

    from pyRGBA.cameras.camera import Camera

    cam = Camera()
    
    cam.zoom(10.0)
    cam.rotateX(np.pi)
    cam.rotateY(np.pi)
    cam.rotateZ(np.pi)
    cam.translateX(-0.5)
    cam.translateY(-0.5)
    cam.translateZ(-0.5)

    view  = cam.get_matrix()

    """

      def __init__(self, matrix = None, surface = None,
                   volume = None, transfer = None,
                  size = (864, 480)) :

            # kernel module
            self.kernel_module = cdh.source_file(
            self.base_directory(__file__, "front_to_back.cu"))

            #kernel function definitions
            self.cast_rays_identifier = self.kernel_module.get_function("front_to_back")
            self.transfer_identifier  = self.kernel_module.get_texref("transfer")
            self.volume_identifier    = self.kernel_module.get_texref("volume")

            self.set_matrix(matrix)

            if (surface == None) :
                  self.set_surface(ArraySurface(size = size))
            else :
                  self.set_surface(surface)


            self.volume = None
            self.set_transfer(transfer)

            self.set_sample_size()
            self.set_max_samples()
            self.set_density_scale()
            self.set_brightness()

      def __call__(self):
            self.cast()

      def get_surface(self) :
          return self.surface.get_array()

      def get_sample_size(self) :
              return self.sample_size

      def get_max_samples(self) :
              return self.max_samples

      def get_density_scale(self) :
              return self.density_scale

      def get_brightness(self):
              return self.brightness

      def set_sample_size(self, size = 0.01) :
              self.sample_size = size

      def set_max_samples(self, max = 5000) :
              self.max_samples = max

      def set_density_scale(self, scale = 0.05) :
              self.density_scale = scale

      def set_brightness(self, brightness = 1.0):
              self.brightness = brightness

      def cast(self):
            w, h = self.surface.bounds

            self.cast_rays_identifier(np.uint64(self.surface.device_ptr),
                                          drv.In(self.matrix),
                                          np.int32(w), np.int32(h),
                                          np.int32(self.max_samples),
                                          np.float32(self.density_scale),
                                          np.float32(self.transfer.offset),
                                          np.float32(self.transfer.scale),
                                          np.float32(self.brightness),
                                          np.float32(self.sample_size),
                                          block=self.block,
                                          grid=self.grid
                                          )
	   

      def set_matrix(self, matrix):
		self.matrix = matrix

      def set_surface(self, surface = None, block_size = 32):
		if surface == None:
			self.surface = None
			return

		self.surface = surface
		self.grid = cdh.block_estimate(block_size, self.surface.bounds)
		self.block = (block_size, block_size, 1)

      def send_volume_to_gpu(self, volume = None) :
            if (volume != None) :
                self.set_volume(Unigrid(array = volume, allocate = True))

      def set_volume(self, volume):
            if volume == None:
                  self.volume = None
                  return

            self.volume = volume
            self.volume_identifier.set_flags(self.volume.flags)
            self.volume_identifier.set_filter_mode(self.volume.filter_mode)
            self.volume_identifier.set_address_mode(0, self.volume.address_mode)
            self.volume_identifier.set_address_mode(1, self.volume.address_mode)
            self.volume_identifier.set_array(self.volume.cuda_volume_array)

      def set_transfer(self, transfer):
		if transfer == None:
			self.transfer = None
			return

		self.transfer = transfer
		self.transfer_identifier.set_flags(self.transfer.flags)
		self.transfer_identifier.set_filter_mode(self.transfer.filter_mode)
		self.transfer_identifier.set_address_mode(0, self.transfer.address_mode)
		self.transfer_identifier.set_array(self.transfer.cuda_transfer_array)

      def base_directory(self, dir, file):
         	base = os.path.dirname(dir)
         	src = os.path.join(base, file)
                return src


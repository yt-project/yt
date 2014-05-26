#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import pycuda.driver as drv
import yt.visualization.volume_rendering.theia.helpers.cuda as cdh
import numpy as np

from yt.visualization.volume_rendering.theia.transfer.linear_transfer import LinearTransferFunction
from yt.visualization.volume_rendering.theia.transfer.helper import *


from yt.visualization.volume_rendering.theia.volumes.unigrid import Unigrid
from yt.visualization.volume_rendering.theia.surfaces.array_surface import ArraySurface
import os

class FrontToBackRaycaster:
      r"""
      This is the python wrapper for the CUDA raycaster. This 
      raycaster can act on unigrid data only. Use yt to map 
      AMR data to unigrid prior to using this algorithm. 


      Parameters
      ----------
      matrix : 4x4 numpy.matriix
          ModelViewMatrix defines view onto volume

      surface: yt.visualization.volume_rendering.theia.surface.array_surface.ArraySurface
          The surface to contain the image information from
          the results of the raycasting

      volume: 3D numpy array Float32
          Scalar values to be volume rendered

      tf : yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction
          Transfer function defining how the raycasting results are 
          to be colored.

      size : Tuple of 2 Floats
         If no surface is specified creates an ArraySurface of the specified size.
 
      """

      def __init__(self, matrix = None, surface = None,
                   volume = None, tf = None,
                  size = (864, 480)) :

            # kernel module
            self.kernel_module = cdh.source_file(
            self.base_directory(__file__, "front_to_back.cu"))

            #kernel function definitions
            self.cast_rays_identifier = \
                self.kernel_module.get_function("front_to_back")
            self.transfer_identifier  = \
                self.kernel_module.get_texref("transfer")
            self.volume_identifier    = \
                self.kernel_module.get_texref("volume")

            self.set_matrix(matrix)

            if (surface == None) :
                  self.set_surface(ArraySurface(size = size))
            else :
                  self.set_surface(surface)


            self.volume = None
            self.set_transfer(tf)

            self.set_sample_size()
            self.set_max_samples()
            self.set_opacity()
            self.set_brightness()

      """
          Envoke the ray caster to cast rays
      """
      def __call__(self):
            self.cast()
      """
          Returns
          ----------
          A  2d numpy array containing the image of the
          volumetric rendering
      """
      def get_surface(self) :
          return self.surface.get_array()

      """
          Returns
          ----------
          The sample size the ray caster is
          configured to take.
      """
      def get_sample_size(self) :
              return self.sample_size

      """
          Returns
          ----------
          The Max number of samples per ray the ray caster is
          configured to take.
      """
      def get_max_samples(self) :
              return self.max_samples

      """
          Returns
          ----------
          The Global density scalar.
      """
      def get_opacity(self) :
              return self.opacity

      """
          Returns
          ----------
          The Global brightness scalar.
         
      """
      def get_brightness(self):
              return self.brightness

      """
          Parameters
          ----------
          size : scalra Float
              The distance between each sample, a smaller size will result in
              more samples and can cause performance loss.
      """
      def set_sample_size(self, size = 0.01) :
              self.sample_size = size

      """
          Parameters
          ----------
          max : scalar Int
              The limit on the number of samples the ray caster will
              take per ray.
      """
      def set_max_samples(self, max = 5000) :
              self.max_samples = max

      """
          Parameters
          ----------
          scale : scalar Float
              Global multiplier on volume data
      """
      def set_opacity(self, scale = 0.05) :
              self.opacity = scale

      """
          Parameters
          ----------
          brightness : scalar Float
              Global multiplier on returned alpha values
      """
      def set_brightness(self, brightness = 1.0):
              self.brightness = brightness

      """
          Causes the ray caster to act on the volume data
          and puts the results in the surface array.
      """
      def cast(self):
            w, h = self.surface.bounds

            self.cast_rays_identifier(np.uint64(self.surface.device_ptr),
                                          drv.In(self.matrix),
                                          np.int32(w), np.int32(h),
                                          np.int32(self.max_samples),
                                          np.float32(self.opacity),
                                          np.float32(self.transfer.offset),
                                          np.float32(self.transfer.scale),
                                          np.float32(self.brightness),
                                          np.float32(self.sample_size),
                                          block=self.block,
                                          grid=self.grid
                                          )
	   

      """
          Set the ModelViewMatrix, this does not send data to the GPU
          Parameters
          ----------
          matrix : 4x4 numpy.matrix
             ModelViewMatrix 
              
      """
      def set_matrix(self, matrix):
		self.matrix = matrix

      """
          Setup the image array for ray casting results
          Parameters
          ----------
          surface : yt.visualization..volume_rendering.theia.surfaces.ArraySurface
              Surface to contain results of the ray casting
      """
      def set_surface(self, surface = None, block_size = 32):
		if surface == None:
			self.surface = None
			return

		self.surface = surface
		self.grid = cdh.block_estimate(block_size, self.surface.bounds)
		self.block = (block_size, block_size, 1)

      """
          This function will convert a 3d numpy array into a unigrid volume
          and move the data to the gpu texture memory
          Parameters
          ----------
          volume : 3D numpy array float32
              Contains scalar volumes to be acted on by the ray caster
          
      """
      def send_volume_to_gpu(self, volume = None) :
            if (volume != None) :
                self.set_volume(Unigrid(array = volume, allocate = True))

      """
          Parameters
          ----------
          volume : yt.visualization..volume_rendering.theia.volumes.Unigrid
              Contains scalar volumes to be acted on by the ray caster
           
      """
      def set_volume(self, volume = None):
            if volume == None:
                  self.volume = None
                  return

            self.volume = volume
            self.volume_identifier.set_flags(self.volume.flags)
            self.volume_identifier.set_filter_mode(self.volume.filter_mode)
            self.volume_identifier.set_address_mode(0, self.volume.address_mode)
            self.volume_identifier.set_address_mode(1, self.volume.address_mode)
            self.volume_identifier.set_array(self.volume.cuda_volume_array)

      """
          Parameters
          ----------
          tf : yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction
             Used to color the results of the raycasting

      """
      def set_transfer(self, tf = None):
		if tf == None:
			self.transfer = None
			return

		self.transfer = LinearTransferFunction(arrays = yt_to_rgba_transfer(tf), range = tf.x_bounds)
		self.transfer_identifier.set_flags(self.transfer.flags)
		self.transfer_identifier.set_filter_mode(self.transfer.filter_mode)
		self.transfer_identifier.set_address_mode(0, self.transfer.address_mode)
		self.transfer_identifier.set_array(self.transfer.cuda_transfer_array)

      """
          Attach the base directory path to the desired source file.
          Parameters
          ----------
          dir : string  
                Directory where source file is located

          file : string 
                Source file name

          
      """
      def base_directory(self, dir, file):
         	base = os.path.dirname(dir)
         	src = os.path.join(base, file)
                return src


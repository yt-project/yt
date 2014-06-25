#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import pycuda.driver   as drv
from   pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray
import pycuda.autoinit
import math
import os

def sample_path():
	#return '/pfs/sw/cuda-5.5/samples/common/inc'
	return os.environ.get('CUDA_SAMPLES')

def cuda_helpers_path():
	return os.path.dirname(__file__) + "/../"

def source_file(path):
	return SourceModule(open(path).read(),
                  include_dirs=[sample_path(),cuda_helpers_path()],
                  no_extern_c=True, keep=True)

def block_estimate(thread_size, shape):
      (x, y) = shape
      return (int(math.ceil(float(x) / thread_size)), int(math.ceil(float(y) / thread_size)), 1)

def np3d_to_device_array(np_array, allow_surface_bind=True):
      d, h, w = np_array.shape

      descr = drv.ArrayDescriptor3D()
      descr.width = w
      descr.height = h
      descr.depth = d
      descr.format = drv.dtype_to_array_format(np_array.dtype)
      descr.num_channels = 1
      descr.flags = 0

      if allow_surface_bind:
            descr.flags = drv.array3d_flags.SURFACE_LDST

      device_array = drv.Array(descr)

      copy = drv.Memcpy3D()
      copy.set_src_host(np_array)
      copy.set_dst_array(device_array)
      copy.width_in_bytes = copy.src_pitch = np_array.strides[1]
      copy.src_height = copy.height = h
      copy.depth = d

      copy()

      return device_array

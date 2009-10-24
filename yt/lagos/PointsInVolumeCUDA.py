from yt.funcs import *
from yt.arraytypes import *
import math, time

import pycuda.driver as cuda
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
cuda.init()

CUDA_BLOCK_SIZE = 22

class VolumeFinder(object):
    def __init__(self, points, res):
        #print cuda.mem_get_info()
        self.res = res
        if not cuda.Device.count() >= 1:
            raise NoCUDAException
        #raise NoCUDAException

        npoints = points.shape[0]

        self.points = []
        for i in range(3):
            pp = points[:, i].astype('float32')
            self.points.append(gpuarray.to_gpu(pp))

        self.valid = gpuarray.zeros((npoints,), 'float32')

        self._mod = cuda.SourceModule(routine)#, keep=True)
        self._func = self._mod.get_function("InVolume")

    # We don't need a delete function, becaue refcounting on the arrays
    # should handle deallocation on the GPU

    def __call__(self, grid):
        t1 = time.time()
        #print cuda.mem_get_info()
        mask = grid.child_mask.astype('float32').copy().ravel()
        ni, nj, nk = grid.ActiveDimensions.astype("int32")
        
        # Now we copy over our new grid mask
        #import pdb;pdb.set_trace()
        mask_gpu = cuda.mem_alloc(mask.size * 4)
        cuda.memcpy_htod(mask_gpu, mask)

        # We need to do left and right, too
        LE, RE = grid.LeftEdge.astype("float32"), grid.RightEdge.astype("float32")

        v1 = self.valid.get()

        gsize = int(math.ceil(float(self.res)/CUDA_BLOCK_SIZE))

        #import pdb; pdb.set_trace()
        self._func( self.points[0].gpudata, self.points[1].gpudata,
                        self.points[2].gpudata,
                    self.valid.gpudata, mask_gpu,
                    LE[0], LE[1], LE[2],
                    RE[0], RE[1], RE[2],
                        just_one(grid['x'].astype('float32')), ni, nj, nk,
                    block=(CUDA_BLOCK_SIZE,CUDA_BLOCK_SIZE,1),
                    grid=(gsize, gsize), time_kernel=True)
        v2 = self.valid.get()

        return na.where(na.logical_xor(v1, v2))[0]



routine = """
__global__ void InVolume(float *points_x, float *points_y, float *points_z,
                         float *valid, float *grid_mask,
                         float left_x, float left_y, float left_z,
                         float right_x, float right_y, float right_z,
                         float dx, int ni, int nj, int nk)
{

    /* Our methodology should be to iterate over every *point* and then mark
       the array with the result of that. We will not do any skipping of
       points, I think? */

  /* My index in the array */
  int pindex = blockIdx.x * blockDim.x + threadIdx.x;

  float point_x = points_x[pindex];
  float point_y = points_y[pindex];
  float point_z = points_z[pindex];

  int addition = (valid[pindex] == 0);
  addition *= (int) ((point_x > left_x) && (point_x < right_x));
  addition *= (int) ((point_y > left_y) && (point_y < right_y));
  addition *= (int) ((point_z > left_z) && (point_z < right_z));
  
  int i, j, k;

  i = (int) ((point_x - left_x) / dx);
  i = min(max(i, 0), ni - 1);

  j = (int) ((point_y - left_y) / dx);
  j = min(max(j, 0), nj - 1);

  k = (int) ((point_z - left_z) / dx);
  k = min(max(k, 0), nk - 1);

  int index = ( ( ( k * nj ) + j ) * ni + i );

  addition *= grid_mask[index];

  valid[pindex] += addition;

    /* Move the trilinear interpolation here, too? */

}
"""

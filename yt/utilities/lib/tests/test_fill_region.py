import numpy as np

from yt.testing import assert_equal
from yt.utilities.lib.misc_utilities import fill_region

NDIM = 32

def test_fill_region():
    for level in range(2):
        rf = 2**level
        output_fields = [np.zeros((NDIM*rf,NDIM*rf,NDIM*rf), "float64")
                         for i in range(3)]
        input_fields = [np.empty(NDIM**3, "float64")
                         for i in range(3)]
        v = np.mgrid[0.0:1.0:NDIM*1j, 0.0:1.0:NDIM*1j, 0.0:1.0:NDIM*1j]
        input_fields[0][:] = v[0].ravel()
        input_fields[1][:] = v[1].ravel()
        input_fields[2][:] = v[2].ravel()
        left_index = np.zeros(3, "int64")
        ipos = np.empty((NDIM**3, 3), dtype="int64")
        ind = np.indices((NDIM,NDIM,NDIM))
        ipos[:,0] = ind[0].ravel()
        ipos[:,1] = ind[1].ravel()
        ipos[:,2] = ind[2].ravel()
        ires = np.zeros(NDIM*NDIM*NDIM, "int64")
        ddims = np.array([NDIM, NDIM, NDIM], dtype="int64") * rf
        fill_region(input_fields, output_fields, level,
                    left_index, ipos, ires, ddims, 2)
        for r in range(level + 1):
            for o, i in zip(output_fields, v):
                assert_equal( o[r::rf,r::rf,r::rf], i)


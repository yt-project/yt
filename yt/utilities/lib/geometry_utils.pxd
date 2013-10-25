"""
Particle Deposition onto Octs




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
cimport numpy as np
cimport cython

DEF ORDER_MAX=20

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t spread_bits(np.uint64_t x):
    # This magic comes from http://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints
    x=(x|(x<<20))&(<np.uint64_t>0x000001FFC00003FF)
    x=(x|(x<<10))&(<np.uint64_t>0x0007E007C00F801F)
    x=(x|(x<<4))&(<np.uint64_t>0x00786070C0E181C3)
    x=(x|(x<<2))&(<np.uint64_t>0x0199219243248649)
    x=(x|(x<<2))&(<np.uint64_t>0x0649249249249249)
    x=(x|(x<<2))&(<np.uint64_t>0x1249249249249249)
    return x

@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *DRE):
    cdef int i
    cdef np.float64_t dds[3]
    cdef np.uint64_t ii[3], mi
    for i in range(3):
        dds[i] = (DRE[i] - DLE[i]) / (1 << ORDER_MAX)
    ii[0] = <np.uint64_t> ((x - DLE[0])/dds[0])
    ii[1] = <np.uint64_t> ((y - DLE[1])/dds[1])
    ii[2] = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = 0
    mi |= spread_bits(ii[2])<<0
    mi |= spread_bits(ii[1])<<1
    mi |= spread_bits(ii[0])<<2
    return mi

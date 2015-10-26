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
DEF INDEX_MAX=2097151

#-----------------------------------------------------------------------------
# 21 bits spread over 64 with 2 bits in between
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t spread_64bits_by2(np.uint64_t x):
    # This magic comes from http://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints
    # Only reversible up to 2097151
    # Select highest 21 bits (Required to be reversible to 21st bit)
    # x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---k jihg fedc ba98 7654 3210
    x=(x&(<np.uint64_t>0x00000000001FFFFF))
    # x = ---- ---- ---- ---- ---- ---k jihg fedc ba-- ---- ---- ---- ---- --98 7654 3210
    x=(x|(x<<20))&(<np.uint64_t>0x000001FFC00003FF)
    # x = ---- ---- ---- -kji hgf- ---- ---- -edc ba-- ---- ---- 9876 5--- ---- ---4 3210
    x=(x|(x<<10))&(<np.uint64_t>0x0007E007C00F801F)
    # x = ---- ---- -kji h--- -gf- ---- -edc ---- ba-- ---- 987- ---6 5--- ---4 32-- --10
    x=(x|(x<<4))&(<np.uint64_t>0x00786070C0E181C3)
    # x = ---- ---k ji-- h--g --f- ---e d--c --b- -a-- --98 --7- -6-- 5--- -43- -2-- 1--0
    x=(x|(x<<2))&(<np.uint64_t>0x0199219243248649)
    # x = ---- -kj- -i-- h--g --f- -e-- d--c --b- -a-- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x=(x|(x<<2))&(<np.uint64_t>0x0649249249249249)
    # x = ---k --j- -i-- h--g --f- -e-- d--c --b- -a-- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x=(x|(x<<2))&(<np.uint64_t>0x1249249249249249)
    return x

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t compact_64bits_by2(np.uint64_t x):
    # Reversed magic
    x=x&(<np.uint64_t>0x1249249249249249)
    x=(x|(x>>2))&(<np.uint64_t>0x0649249249249249)
    x=(x|(x>>2))&(<np.uint64_t>0x0199219243248649)
    x=(x|(x>>2))&(<np.uint64_t>0x00786070C0E181C3)
    x=(x|(x>>4))&(<np.uint64_t>0x0007E007C00F801F)
    x=(x|(x>>10))&(<np.uint64_t>0x000001FFC00003FF)
    x=(x|(x>>20))&(<np.uint64_t>0x00000000001FFFFF) 
    return x

#-----------------------------------------------------------------------------
# 10 bits spread over 32 with 2 bits in between
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint32_t spread_32bits_by2(np.uint32_t x):
    # Only reversible up to 1023
    # Select highest 10 bits (Required to be reversible to 10st bit)
    # x = ---- ---- ---- ---- ---- --98 7654 3210
    x=(x&(<np.uint32_t>0x000003FF))
    # x = ---- --98 ---- ---- ---- ---- 7654 3210
    x=(x|(x<<16))&(<np.uint32_t>0xFF0000FF)
    # x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x=(x|(x<<8))&(<np.uint32_t>0x0300F00F)
    # x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x=(x|(x<<4))&(<np.uint32_t>0x030C30C3)
    # x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x=(x|(x<<2))&(<np.uint32_t>0x09249249)
    return x

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint32_t compact_32bits_by2(np.uint32_t x):
    # Reversed magic
    x=x&(<np.uint32_t>0x09249249)
    x=(x|(x>>2))&(<np.uint32_t>0x030C30C3)
    x=(x|(x>>4))&(<np.uint32_t>0x0300F00F)
    x=(x|(x>>8))&(<np.uint32_t>0xFF0000FF)
    x=(x|(x>>16))&(<np.uint32_t>0x000003FF)
    return x

#-----------------------------------------------------------------------------
# Shortcuts for default
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t spread_bits(np.uint64_t x):
    return spread_64bits_by2(x)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t compact_bits(np.uint64_t x):
    return compact_64bits_by2(x)

@cython.cdivision(True)
cdef inline np.uint64_t encode_morton_64bit(np.uint64_t x_ind, np.uint64_t y_ind, np.uint64_t z_ind):
    cdef np.uint64_t mi
    mi = 0
    mi |= spread_64bits_by2(z_ind)<<0
    mi |= spread_64bits_by2(y_ind)<<1
    mi |= spread_64bits_by2(x_ind)<<2
    return mi

# @cython.cdivision(True)
# cdef inline np.uint64_t decode_morton_64bit(np.uint64_t mi):
#     cdef np.uint64_t ii[3]
#     ii[0] = compact_64bits_by2(mi>>2)
#     ii[1] = compact_64bits_by2(mi>>1)
#     ii[2] = compact_64bits_by2(mi>>0)
#     return ii

@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *DRE, np.int32_t order):
    cdef int i
    cdef np.float64_t dds[3]
    cdef np.uint64_t ii[3]
    cdef np.uint64_t mi
    for i in range(3):
        dds[i] = (DRE[i] - DLE[i]) / (1 << order)
    ii[0] = <np.uint64_t> ((x - DLE[0])/dds[0])
    ii[1] = <np.uint64_t> ((y - DLE[1])/dds[1])
    ii[2] = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = 0
    mi |= spread_bits(ii[2])<<0
    mi |= spread_bits(ii[1])<<1
    mi |= spread_bits(ii[0])<<2
    return mi
  
# This dosn't seem to be much, if at all, faster...
@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_mml(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *dds):
    cdef int i
    cdef np.uint64_t ii[3]
    cdef np.uint64_t mi
    ii[0] = <np.uint64_t> ((x - DLE[0])/dds[0])
    ii[1] = <np.uint64_t> ((y - DLE[1])/dds[1])
    ii[2] = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = 0
    mi |= spread_bits(ii[2])<<0
    mi |= spread_bits(ii[1])<<1
    mi |= spread_bits(ii[0])<<2
    return mi

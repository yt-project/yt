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
from libc.float cimport DBL_MANT_DIG
from libc.math cimport frexp,ldexp

DEF ORDER_MAX=20
DEF INDEX_MAX=2097151
    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.int64_t ifrexp(np.float64_t x, np.int64_t *e):
    cdef np.float64_t m
    cdef int e0
    m = frexp(x,&e0)
    e[0] = <np.int64_t>e0
    return <np.int64_t>ldexp(m,DBL_MANT_DIG)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.int64_t msdb(np.int64_t a, np.int64_t b):
    """Get the most significant differing bit between a and b."""
    cdef np.int64_t c, ndx
    c = a ^ b
    ndx = 0
    while (1 < c):
        c = (c >> 1)
        ndx+=1
    return ndx

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.int64_t xor_msb(np.float64_t a, np.float64_t b):
    """Get the exponent of the highest differing bit between a and b"""
    # Get mantissa and exponents for each number
    cdef np.int64_t a_m, a_e, b_m, b_e, x, y, z
    b_e = 0
    a_e = 0
    a_m = ifrexp(a,&a_e)
    b_m = ifrexp(b,&b_e)
    x = a_e
    y = b_e
    # Compare mantissa if exponents equal
    if x == y:
        z = msdb(a_m,b_m)
        x = x - z
        return x
    # Otherwise return largest exponent
    if y < x:
        return x
    else:
        return y

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int compare_floats_morton(np.float64_t p[3], np.float64_t q[3]):
    cdef int j, out, dim
    cdef np.int64_t x, y
    x = -308
    y = 0
    dim = 0
    for j in range(3):#[::-1]:
        y = xor_msb(p[j],q[j])
        if x < y:
           x = y
           dim = j
    if p[dim] < q[dim]:
        out = 1
    else:
        out = 0
    return out

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.float64_t euclidean_distance(np.float64_t p[3], np.float64_t q[3]):
    cdef int j
    cdef np.float64_t d
    d = 0.0
    for j in range(3):
        d+=(p[j]-q[j])**2
    return np.sqrt(d)

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
# Shortcuts for default (May be invalid to call inline from within inline...)
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

@cython.cdivision(True)
cdef inline void decode_morton_64bit(np.uint64_t mi, np.uint64_t *p):
    p[0] = compact_64bits_by2(mi>>2)
    p[1] = compact_64bits_by2(mi>>1)
    p[2] = compact_64bits_by2(mi>>0)

@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *DRE, np.int32_t order):
    cdef int i
    cdef np.float64_t dds[3]
    cdef np.uint64_t x_ind, y_ind, z_ind
    cdef np.uint64_t mi
    for i in range(3):
        dds[i] = (DRE[i] - DLE[i]) / (1 << order)
    x_ind = <np.uint64_t> ((x - DLE[0])/dds[0])
    y_ind = <np.uint64_t> ((y - DLE[1])/dds[1])
    z_ind = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = encode_morton_64bit(x_ind,y_ind,z_ind)
    return mi
  
# This dosn't seem to be much, if at all, faster...
@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_mml(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *dds):
    cdef int i
    cdef np.uint64_t x_ind, y_ind, z_ind
    cdef np.uint64_t mi
    x_ind = <np.uint64_t> ((x - DLE[0])/dds[0])
    y_ind = <np.uint64_t> ((y - DLE[1])/dds[1])
    z_ind = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = encode_morton_64bit(x_ind,y_ind,z_ind)
    return mi

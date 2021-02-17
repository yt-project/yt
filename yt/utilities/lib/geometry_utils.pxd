"""
Particle Deposition onto Octs




"""

cimport cython
cimport numpy as np
from libc.float cimport DBL_MANT_DIG
from libc.math cimport frexp, ldexp, sqrt

DEF ORDER_MAX=20
DEF INDEX_MAX_64=2097151
# TODO: Handle error for indices past max
DEF XSHIFT=2
DEF YSHIFT=1
DEF ZSHIFT=0

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.int64_t ifrexp(np.float64_t x, np.int64_t *e):
    cdef np.float64_t m
    cdef int e0 = 0
    m = frexp(x,&e0)
    e[0] = <np.int64_t>e0
    return <np.int64_t>ldexp(m,<int>DBL_MANT_DIG)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.int64_t msdb(np.int64_t a, np.int64_t b):
    """Get the most significant differing bit between a and b."""
    cdef np.int64_t c, ndx
    c = a ^ b
    ndx = 0
    while (0 < c):
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
    x = <np.int64_t> ((a_e+1)*DBL_MANT_DIG)
    y = <np.int64_t> ((b_e+1)*DBL_MANT_DIG)
    # Compare mantissa if exponents equal
    if x == y:
        if a_m == b_m: return 0
        z = msdb(a_m,b_m)
        #if 1: return z
        x = x - z
        return x-1 # required so that xor_msb(0.0,1.0)!=xor_msb(1.0,1.0)
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
    x = -9999999999
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
cdef inline np.float64_t euclidean_distance(np.float64_t[:] p, np.float64_t[:] q):
    cdef int j
    cdef np.float64_t d
    d = 0.0
    for j in range(3):
        d+=(p[j]-q[j])**2
    return sqrt(d)

# Todo: allow radius reported independently in each dimension for rectangular domain
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.float64_t smallest_quadtree_box(np.float64_t p[3], np.float64_t q[3], np.int32_t order,
                                               np.float64_t DLE[3], np.float64_t DRE[3],
                                               np.float64_t *cx, np.float64_t *cy, np.float64_t *cz):
    cdef int j
    cdef np.float64_t c[3]
    cdef np.uint64_t pidx[3]
    # cdef np.uint64_t qidx[3]
    for j in range(3):
        pidx[j] = 0
        # qidx[j] = 0
    cdef np.uint64_t pidx_next[3]
    cdef np.uint64_t qidx_next[3]
    cdef np.float64_t dds[3]
    cdef np.float64_t rad
    cdef int lvl = 0
    cdef int done = 0
    while not done:
        if (lvl+1 >= order):
            done = 1
        for j in range(3):
            dds[j] = (DRE[j] - DLE[j])/(1 << (<int> lvl+1))
            pidx_next[j] = <np.uint64_t>((p[j] - DLE[j])/dds[j])
            qidx_next[j] = <np.uint64_t>((q[j] - DLE[j])/dds[j])
        for j in range(3):
            if pidx_next[j]!=qidx_next[j]:
                done = 1
                break
        if not done:
            for j in range(3):
                pidx[j] = pidx_next[j]
                # qidx[j] = qidx_next[j]
            lvl+=1
    rad = 0.0
    for j in range(3):
        dds[j] = (DRE[j] - DLE[j])/(1 << lvl)
        c[j] = dds[j]*(<np.float64_t>pidx[j]+0.5)
        rad+=((dds[j]/2.0)**2)
    cx[0] = c[0]
    cy[0] = c[1]
    cz[0] = c[2]
    return sqrt(rad)

#-----------------------------------------------------------------------------
# 21 bits spread over 64 with 3 bits in between
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t spread_64bits_by3(np.uint64_t x):
    x=(x&(<np.uint64_t>0x00000000001FFFFF))
    x=(x|(x<<20))*(<np.uint64_t>0x000001FFC00003FF)

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

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t masked_merge_64bit(np.uint64_t a, np.uint64_t b, np.uint64_t mask):
    # https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
    return a ^ ((a ^ b) & mask)

@cython.cdivision(True)
cdef inline np.uint64_t encode_morton_64bit(np.uint64_t x_ind, np.uint64_t y_ind, np.uint64_t z_ind):
    cdef np.uint64_t mi
    mi = 0
    mi |= spread_64bits_by2(z_ind)<<ZSHIFT
    mi |= spread_64bits_by2(y_ind)<<YSHIFT
    mi |= spread_64bits_by2(x_ind)<<XSHIFT
    return mi

@cython.cdivision(True)
cdef inline void decode_morton_64bit(np.uint64_t mi, np.uint64_t *p):
    p[0] = compact_64bits_by2(mi>>XSHIFT)
    p[1] = compact_64bits_by2(mi>>YSHIFT)
    p[2] = compact_64bits_by2(mi>>ZSHIFT)

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

@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_relative(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *DRE,
                               np.int32_t order1, np.int32_t order2):
    cdef int i
    cdef np.float64_t dds1[3]
    cdef np.float64_t dds2[3]
    cdef np.float64_t DLE2[3]
    cdef np.uint64_t x_ind, y_ind, z_ind
    cdef np.uint64_t mi2
    for i in range(3):
        dds1[i] = (DRE[i] - DLE[i]) / (1 << order1)
        dds2[i] = dds1[i] / (1 << order2)
    DLE2[0] = <np.float64_t> (<np.uint64_t> ((x - DLE[0])/dds1[0])) * dds1[0]
    DLE2[1] = <np.float64_t> (<np.uint64_t> ((y - DLE[1])/dds1[1])) * dds1[1]
    DLE2[2] = <np.float64_t> (<np.uint64_t> ((z - DLE[2])/dds1[2])) * dds1[2]
    x_ind = <np.uint64_t> ((x - DLE2[0])/dds2[0])
    y_ind = <np.uint64_t> ((y - DLE2[1])/dds2[1])
    z_ind = <np.uint64_t> ((z - DLE2[2])/dds2[2])
    mi2 = encode_morton_64bit(x_ind,y_ind,z_ind)
    return mi2


# This dosn't seem to be much, if at all, faster...
@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_dds(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *dds):
    cdef np.uint64_t x_ind, y_ind, z_ind
    cdef np.uint64_t mi
    x_ind = <np.uint64_t> ((x - DLE[0])/dds[0])
    y_ind = <np.uint64_t> ((y - DLE[1])/dds[1])
    z_ind = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = encode_morton_64bit(x_ind,y_ind,z_ind)
    return mi

@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_relative_dds(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *dds1, np.float64_t *dds2):
    cdef np.float64_t DLE2[3]
    cdef np.uint64_t x_ind, y_ind, z_ind
    cdef np.uint64_t mi2
    DLE2[0] = <np.float64_t> (<np.uint64_t> ((x - DLE[0])/dds1[0])) * dds1[0]
    DLE2[1] = <np.float64_t> (<np.uint64_t> ((y - DLE[1])/dds1[1])) * dds1[1]
    DLE2[2] = <np.float64_t> (<np.uint64_t> ((z - DLE[2])/dds1[2])) * dds1[2]
    x_ind = <np.uint64_t> ((x - DLE2[0])/dds2[0])
    y_ind = <np.uint64_t> ((y - DLE2[1])/dds2[1])
    z_ind = <np.uint64_t> ((z - DLE2[2])/dds2[2])
    mi2 = encode_morton_64bit(x_ind,y_ind,z_ind)
    return mi2


@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_split_dds(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *dds, np.uint64_t *p):
    cdef np.uint64_t mi
    p[0] = <np.uint64_t> ((x - DLE[0])/dds[0])
    p[1] = <np.uint64_t> ((y - DLE[1])/dds[1])
    p[2] = <np.uint64_t> ((z - DLE[2])/dds[2])
    mi = encode_morton_64bit(p[0], p[1], p[2])
    return mi

@cython.cdivision(True)
cdef inline np.uint64_t bounded_morton_split_relative_dds(np.float64_t x, np.float64_t y, np.float64_t z,
                               np.float64_t *DLE, np.float64_t *dds1, np.float64_t *dds2,
                               np.uint64_t *p2):
    cdef np.float64_t DLE2[3]
    cdef np.uint64_t mi2
    DLE2[0] = DLE[0] + <np.float64_t> (<np.uint64_t> ((x - DLE[0])/dds1[0])) * dds1[0]
    DLE2[1] = DLE[1] + <np.float64_t> (<np.uint64_t> ((y - DLE[1])/dds1[1])) * dds1[1]
    DLE2[2] = DLE[2] + <np.float64_t> (<np.uint64_t> ((z - DLE[2])/dds1[2])) * dds1[2]
    p2[0] = <np.uint64_t> ((x - DLE2[0])/dds2[0])
    p2[1] = <np.uint64_t> ((y - DLE2[1])/dds2[1])
    p2[2] = <np.uint64_t> ((z - DLE2[2])/dds2[2])
    mi2 = encode_morton_64bit(p2[0], p2[1], p2[2])
    return mi2


cdef np.uint32_t morton_neighbors_coarse(np.uint64_t mi1, np.uint64_t max_index1,
                                         bint periodicity[3], np.uint32_t nn,
                                         np.uint32_t[:,:] index,
                                         np.uint64_t[:,:] ind1_n,
                                         np.uint64_t[:] neighbors)

cdef np.uint32_t morton_neighbors_refined(np.uint64_t mi1, np.uint64_t mi2,
                                          np.uint64_t max_index1, np.uint64_t max_index2,
                                          bint periodicity[3], np.uint32_t nn,
                                          np.uint32_t[:,:] index,
                                          np.uint64_t[:,:] ind1_n,
                                          np.uint64_t[:,:] ind2_n,
                                          np.uint64_t[:] neighbors1,
                                          np.uint64_t[:] neighbors2)

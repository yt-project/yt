"""
Bit array functions



"""


import numpy as np

cimport cython
cimport numpy as np


cdef inline void ba_set_value(np.uint8_t *buf, np.uint64_t ind,
                              np.uint8_t val) noexcept nogil:
    # This assumes 8 bit buffer.  If value is greater than zero (thus allowing
    # us to use 1-255 as 'True') then we identify first the index in the buffer
    # we are setting.  We do this by truncating the index by bit-shifting to
    # the left three times, essentially dividing it by eight (and taking the
    # floor.)
    # The next step is to turn *on* what we're attempting to turn on, which
    # means taking our index and truncating it to the first 3 bits (which we do
    # with an & operation) and then turning on the correct bit.
    #
    # So if we're asking for index 33 in the bitarray, we would want the 4th
    # uint8 element, then the 2nd bit (index 1).
    #
    # To turn it on, we logical *or* with that.  To turn it off, we logical
    # *and* with the *inverse*, which will allow everything *but* that bit to
    # stay on.
    if val > 0:
        buf[ind >> 3] |= (1 << (ind & 7))
    else:
        buf[ind >> 3] &= ~(1 << (ind & 7))

cdef inline np.uint8_t ba_get_value(np.uint8_t *buf, np.uint64_t ind) noexcept nogil:
    cdef np.uint8_t rv = (buf[ind >> 3] & (1 << (ind & 7)))
    if rv == 0: return 0
    return 1

cdef inline void ba_set_range(np.uint8_t *buf, np.uint64_t start_ind,
                              np.uint64_t stop_ind, np.uint8_t val) nogil:
    # Should this be inclusive of both end points?  I think it should not, to
    # match slicing semantics.
    #
    # We need to figure out the first and last values, and then we just set the
    # ones in-between to 255.
    if stop_ind < start_ind: return
    cdef np.uint64_t i
    cdef np.uint8_t j, bitmask
    cdef np.uint64_t buf_start = start_ind >> 3
    cdef np.uint64_t buf_stop = stop_ind >> 3
    cdef np.uint8_t start_j = start_ind & 7
    cdef np.uint8_t stop_j = stop_ind & 7
    if buf_start == buf_stop:
        for j in range(start_j, stop_j):
            ba_set_value(&buf[buf_start], j, val)
        return
    bitmask = 0
    for j in range(start_j, 8):
        bitmask |= (1 << j)
    if val > 0:
        buf[buf_start] |= bitmask
    else:
        buf[buf_start] &= ~bitmask
    if val > 0:
        bitmask = 255
    else:
        bitmask = 0
    for i in range(buf_start + 1, buf_stop):
        buf[i] = bitmask
    bitmask = 0
    for j in range(0, stop_j):
        bitmask |= (1 << j)
    if val > 0:
        buf[buf_stop] |= bitmask
    else:
        buf[buf_stop] &= ~bitmask


cdef inline np.uint8_t _num_set_bits( np.uint8_t b ):
    # https://stackoverflow.com/questions/30688465/how-to-check-the-number-of-set-bits-in-an-8-bit-unsigned-char
    b = b - ((b >> 1) & 0x55)
    b = (b & 0x33) + ((b >> 2) & 0x33)
    return (((b + (b >> 4)) & 0x0F) * 0x01)

cdef class bitarray:
    cdef np.uint8_t *buf
    cdef np.uint64_t size
    cdef np.uint64_t buf_size # Not exactly the same
    cdef np.uint8_t final_bitmask
    cdef public object ibuf

    cdef void _set_value(self, np.uint64_t ind, np.uint8_t val)
    cdef np.uint8_t _query_value(self, np.uint64_t ind)
    cdef void _set_range(self, np.uint64_t start, np.uint64_t stop, np.uint8_t val)
    cdef np.uint64_t _count(self)
    cdef bitarray _logical_and(self, bitarray other, bitarray result = *)
    cdef bitarray _logical_or(self, bitarray other, bitarray result = *)
    cdef bitarray _logical_xor(self, bitarray other, bitarray result = *)

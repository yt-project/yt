"""
Simple readers for fortran unformatted data, specifically for the Tiger code.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
cimport numpy as np
cimport cython

from stdio cimport fopen, fclose, FILE

cdef extern from "stdio.h":
    cdef int SEEK_SET
    cdef int SEEK_CUR
    cdef int SEEK_END
    int fseek(FILE *stream, long offset, int whence)
    size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
    long ftell(FILE *stream)

@cython.boundscheck(False)
@cython.wraparound(False)
def read_tiger_section(
                     char *fn,
                     np.ndarray[np.int64_t, ndim=1] slab_start,
                     np.ndarray[np.int64_t, ndim=1] slab_size,
                     np.ndarray[np.int64_t, ndim=1] root_size,
                     int offset = 36):
    cdef int strides[3]
    strides[0] = 1
    strides[1] = root_size[0] * strides[0]
    strides[2] = strides[1] * root_size[1] + 2
    cdef np.int64_t i, j, k
    cdef np.ndarray buffer = np.zeros(slab_size, dtype='float32', order='F')
    cdef FILE *f = fopen(fn, "rb")
    #for i in range(3): offset += strides[i] * slab_start[i]
    cdef np.int64_t pos = 0
    cdef np.int64_t moff = 0
    cdef float *data = <float *> buffer.data
    fseek(f, offset, 0)
    # If anybody wants to convert this loop to a SEEK_CUR, that'd be great.
    for i in range(slab_size[2]):
        for j in range(slab_size[1]):
            moff = (slab_start[0]    ) * strides[0] \
                 + (slab_start[1] + j) * strides[1] \
                 + (slab_start[2] + i) * strides[2]
            #print offset + 4 * moff, pos
            fseek(f, offset + 4 * moff, SEEK_SET)
            fread(<void *> (data + pos), 4, slab_size[0], f)
            pos += slab_size[0]
    return buffer

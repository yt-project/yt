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

#cdef inline int imax(int i0, int i1):
    #if i0 > i1: return i0
    #return i1

cdef extern from "endian_swap.h":
    void FIX_SHORT( unsigned short )
    void FIX_LONG( unsigned )
    void FIX_FLOAT( float )

cdef extern from "alloca.h":
    void *alloca(int)

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

def count_art_octs(char *fn, long offset,
                   int min_level, int max_level,
                   level_info):
    cdef int nchild = 8
    cdef int i, Lev, next_record, nLevel
    cdef int dummy_records[9]
    cdef int readin
    cdef FILE *f = fopen(fn, "rb")
    fseek(f, offset, SEEK_SET)
    for Lev in range(min_level + 1, max_level + 1):
        fread(dummy_records, sizeof(int), 2, f);
        fread(&nLevel, sizeof(int), 1, f); FIX_LONG(nLevel)
        level_info.append(nLevel)
        fread(dummy_records, sizeof(int), 2, f);
        fread(&next_record, sizeof(int), 1, f); FIX_LONG(next_record)
        print "Record size is:", next_record
        # Offset for one record header we just read
        next_record = (nLevel * (next_record + 2*sizeof(int))) - sizeof(int)
        fseek(f, next_record, SEEK_CUR)
        # Now we skip the second section 
        fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
        next_record = (2*sizeof(int) + readin) * (nLevel * nchild)
        next_record -= sizeof(int)
        fseek(f, next_record, SEEK_CUR)
    fclose(f)

def read_art_tree(char *fn, long offset,
                  int min_level, int max_level,
                  np.ndarray[np.int64_t, ndim=2] oct_indices,
                  np.ndarray[np.int64_t, ndim=1] oct_levels,
                  np.ndarray[np.int64_t, ndim=1] oct_parents,
                  np.ndarray[np.int64_t, ndim=1] oct_mask):
    # This accepts the filename of the ART header and an integer offset that
    # points to the start of the record *following* the reading of iOctFree and
    # nOct.  For those following along at home, we only need to read:
    #   iOctPr, iOctLv
    cdef int nchild = 8
    cdef int i, Lev, cell_ind, iOct, nLevel, nLevCells, ic1, next_record
    cdef int idc, cm
    cdef int iOctPs[3]
    cdef int dummy_records[9]
    cdef int readin
    cdef FILE *f = fopen(fn, "rb")
    fseek(f, offset, SEEK_SET)
    cdef int Level
    cdef int * iNOLL = <int *> alloca(sizeof(int)*(max_level-min_level+1))
    cdef int * iHOLL = <int *> alloca(sizeof(int)*(max_level-min_level+1))
    cell_ind = 0
    cdef int total_cells = 0, total_masked
    cdef int iOctMax = 0
    idc = 0
    for Lev in range(min_level + 1, max_level + 1):
        fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
        fread(&Level, sizeof(int), 1, f); FIX_LONG(Level)
        fread(&iNOLL[Level], sizeof(int), 1, f); FIX_LONG(iNOLL[Level])
        fread(&iHOLL[Level], sizeof(int), 1, f); FIX_LONG(iHOLL[Level])
        fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
        iOct = iHOLL[Level] - 1
        nLevel = iNOLL[Level]
        print "Reading Hierarchy for Level", Lev, Level, nLevel, iOct
        total_cells += nLevel
        for ic1 in range(nLevel):
            iOctMax = imax(iOctMax, iOct)
            #print readin, iOct, nLevel, sizeof(int) 
            next_record = ftell(f)
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            next_record += readin + 2*sizeof(int)
            fread(iOctPs, sizeof(int), 3, f)
            FIX_LONG(iOctPs[0]); FIX_LONG(iOctPs[1]); FIX_LONG(iOctPs[2])
            oct_indices[iOct, 0] = iOctPs[0]
            oct_indices[iOct, 1] = iOctPs[1]
            oct_indices[iOct, 2] = iOctPs[2]
            fread(dummy_records, sizeof(int), 6, f) # skip Nb
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            oct_parents[iOct] = readin - 1
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            oct_levels[iOct] = readin
            fread(&iOct, sizeof(int), 1, f); FIX_LONG(iOct);
            iOct -= 1
            fseek(f, next_record, SEEK_SET)
        total_masked = 0
        for ic1 in range(nLevel * nchild):
            fread(&next_record, sizeof(int), 1, f); FIX_LONG(next_record)
            fread(&idc, sizeof(int), 1, f); FIX_LONG(idc); idc -= 1 + (128**3)
            fread(&cm, sizeof(int), 1, f); FIX_LONG(cm)
            if cm == 0: oct_mask[idc] = 1
            else: total_masked += 1
            fseek(f, next_record - sizeof(int), SEEK_CUR)
        print "Masked cells", total_masked
    fclose(f)
    print "Read this many cells", total_cells, iOctMax

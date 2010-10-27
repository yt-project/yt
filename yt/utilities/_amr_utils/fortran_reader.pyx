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
                   int nhydro_vars,                   
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
        print level_info
        level_info.append(nLevel)
        fread(dummy_records, sizeof(int), 2, f);
        fread(&next_record, sizeof(int), 1, f); FIX_LONG(next_record)
        print "Record size is:", next_record
        # Offset for one record header we just read
        next_record = (nLevel * (next_record + 2*sizeof(int))) - sizeof(int)
        fseek(f, next_record, SEEK_CUR)
        # Now we skip the second section 
        fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
        nhydro_vars = next_record/4-2-3 #nhvar in daniel's code
        #record length is normally 2 pad bytes, 8 + 2 hvars (the 2 is nchem)
        # and then 3 vars, but we can find nhvars only here and not in other
        # file headers 
        next_record = (2*sizeof(int) + readin) * (nLevel * nchild)
        next_record -= sizeof(int)
        fseek(f, next_record, SEEK_CUR)
    print "nhvars",nhydro_vars
    fclose(f)

def read_art_tree(char *fn, long offset,
                  int min_level, int max_level, 
                  np.ndarray[np.int64_t, ndim=2] oct_indices,
                  np.ndarray[np.int64_t, ndim=1] oct_levels,
                  np.ndarray[np.int64_t, ndim=1] oct_parents,
                  np.ndarray[np.int64_t, ndim=1] oct_mask,
                  np.ndarray[np.int64_t, ndim=2] oct_info):
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
    level_offsets = [0]
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
            iOctMax = max(iOctMax, iOct)
            #print readin, iOct, nLevel, sizeof(int) 
            next_record = ftell(f)
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            next_record += readin + 2*sizeof(int)
            fread(iOctPs, sizeof(int), 3, f)
            FIX_LONG(iOctPs[0]); FIX_LONG(iOctPs[1]); FIX_LONG(iOctPs[2])
            oct_indices[iOct, 0] = iOctPs[0]
            oct_indices[iOct, 1] = iOctPs[1]
            oct_indices[iOct, 2] = iOctPs[2]
            oct_info[iOct, 1] = ic1
            #grid_info[iOct, 2] = iOctPr # we don't seem to need this
            fread(dummy_records, sizeof(int), 6, f) # skip Nb
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            oct_parents[iOct] = readin - 1
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            oct_levels[iOct] = readin
            fread(&iOct, sizeof(int), 1, f); FIX_LONG(iOct);
            iOct -= 1
            fseek(f, next_record, SEEK_SET)
        total_masked = 0
        level_offsets.append(ftell(f))
        for ic1 in range(nLevel * nchild):
            fread(&next_record, sizeof(int), 1, f); FIX_LONG(next_record)
            fread(&idc, sizeof(int), 1, f); FIX_LONG(idc); idc -= 1 + (128**3)
            fread(&cm, sizeof(int), 1, f); FIX_LONG(cm)
            if cm == 0: oct_mask[idc] = 1
            else: total_masked += 1
            fseek(f, next_record - sizeof(int), SEEK_CUR)
    fclose(f)
    return level_offsets

def read_art_root_vars(char *fn, long grid_id, long ncell,
                    long root_grid_offset, int nhydro_vars,
                    fields,
                    np.ndarray[np.int64_t, ndim=1] level_info,
                    var):

    cdef FILE *f = fopen(fn, "rb")
    cdef long offset,offset_hvars,offset_var
    cdef int j,l
    cdef float temp
    offset=0
    offset_hvars=0
    offset_var= 0
    #var = np.zeros(len(fields))
    #print 'blah'
    offset += root_grid_offset
    offset += sizeof(int) #pad
    offset += ncell*sizeof(int) #skip iOctCh
    offset += sizeof(int) #pad
    #print 'blah'

    l=0
    offset_hvars += sizeof(int) #pad
    offset_hvars += nhydro_vars*grid_id*sizeof(int) #skip to right before our hvar
    fseek(f,offset+offset_hvars,SEEK_SET)    
    for j in range(nhydro_vars):
        fread(&temp, sizeof(float), 1, f); 
        if j in fields:
            #print 'j',j,'l',l
            FIX_FLOAT(temp)
            #print temp
            #print var
            var[l]=temp
            l+=1
    
    l=0
    offset_var += sizeof(int) #pad
    offset_var += nhydro_vars*ncell*sizeof(int) #skip hvars
    offset_var += sizeof(int) #pad
    offset_var += sizeof(int) #pad
    offset_var += 2*grid_id*sizeof(int) #skip to right before our var
    fseek(f,offset+offset_var,SEEK_SET)    
    for j in range(2):
        fread(&temp, sizeof(float), 1, f); 
        if j+nhydro_vars in fields:
            #print 'blah',j
            FIX_FLOAT(temp)
            var[l]=temp
            l+=1
    fclose(f)

def read_art_vars(char *fn, 
                    int min_level, int max_level, int nhydro_vars, 
                    int grid_level,long grid_id,long child_offset,
                    fields,
                    np.ndarray[np.int64_t, ndim=1] level_offsets,
                    var):
    # nhydro_vars is the number of columns- 3 (adjusting for vars)
    # this is normally 10=(8+2chem species)
    cdef FILE *f = fopen(fn, "rb")
    cdef int record_size = 2+1+1+nhydro_vars+2
    cdef float temp
    cdef int padding[3]
    cdef long offset = 8*grid_id*record_size*sizeof(float)
    fseek(f, level_offsets[grid_level] + offset, SEEK_SET)
    for j in range(8): #iterate over the children
        l = 0
        fread(padding, sizeof(int), 3, f);
        FIX_LONG(padding[0])
        #print "Record Size", padding[0]
        for k in range(nhydro_vars): #iterate over the record
            fread(&temp, sizeof(float), 1, f); FIX_FLOAT(temp)
            #print k, temp
            if k in fields:
                var[j,l] = temp
                l += 1
    fclose(f)

def read_art_grid(int varindex, long ncell,
              np.ndarray[np.int64_t, ndim=1] start_index,
              np.ndarray[np.int32_t, ndim=1] grid_dims,
              np.ndarray[np.float64_t, ndim=3] data,
              np.ndarray[np.int32_t, ndim=3] filled,
              int level, int ref_factor,
              component_grid_info,
              char *fn,                                       
              int min_level, int max_level, int nhydro_vars,
              long root_grid_offset, long child_offset, 
              np.ndarray[np.int64_t, ndim=1] level_offsets):
    cdef int gi, i, j, k, domain, offset, grid_id
    cdef int ir, jr, kr
    cdef int offi, offj, offk, odind
    cdef np.int64_t di, dj, dk
    cdef np.ndarray[np.int64_t, ndim=1] ogrid_info
    cdef np.ndarray[np.int64_t, ndim=1] og_start_index
    cdef np.float64_t temp_data
    cdef np.int64_t end_index[3]
    cdef int to_fill = 0
    # Note that indexing into a cell is:
    #   (k*2 + j)*2 + i
    for i in range(3):
        end_index[i] = start_index[i] + grid_dims[i]
    for gi in range(len(component_grid_info)):
        ogrid_info = component_grid_info[gi]
        domain = ogrid_info[0]
        #print "Loading", domain, ogrid_info
        grid_id = ogrid_info[1]
        og_start_index = ogrid_info[3:]
        for i in range(2*ref_factor):
            di = i + og_start_index[0] * ref_factor
            if di < start_index[0] or di >= end_index[0]: continue
            ir = <int> (i / ref_factor)
            for j in range(2 * ref_factor):
                dj = j + og_start_index[1] * ref_factor
                if dj < start_index[1] or dj >= end_index[1]: continue
                jr = <int> (j / ref_factor)
                for k in range(2 * ref_factor):
                    dk = k + og_start_index[2] * ref_factor
                    if dk < start_index[2] or dk >= end_index[2]: continue
                    kr = <int> (k / ref_factor)
                    offi = di - start_index[0]
                    offj = dj - start_index[1]
                    offk = dk - start_index[2]
                    #print offi, filled.shape[0],
                    #print offj, filled.shape[1],
                    #print offk, filled.shape[2]
                    if filled[offi, offj, offk] == 1: continue

                    odind = (kr*2 + jr)*2 + ir
                    # Replace with an ART-specific reader
                    #temp_data = local_hydro_data.m_var_array[
                    #        level][8*offset + odind]
                    
                    if level > 0:
                        hvars = np.zeros((8,1))
                        read_art_vars(fn, min_level, max_level, nhydro_vars, 
                            level, grid_id, child_offset, [varindex],
                            level_offsets, hvars)
                        temp_data = hvars[odind][0]
                    else:
                        hvars = np.zeros((1))
                        read_art_root_vars(fn, grid_id, ncell, 
                            root_grid_offset, nhydro_vars, [varindex], 
                            level_offsets, hvars)
                        temp_data = hvars[odind]
                    #print temp_data
                    data[offi, offj, offk] = temp_data
                    filled[offi, offj, offk] = 1
                    to_fill += 1
    return to_fill

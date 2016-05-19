"""
Simple readers for fortran unformatted data, specifically for the Tiger code.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

from libc.stdio cimport fopen, fclose, FILE
cimport libc.stdlib as stdlib

#cdef inline int imax(int i0, int i1):
    #if i0 > i1: return i0
    #return i1

cdef extern from "endian_swap.h":
    void FIX_SHORT( unsigned short )
    void FIX_LONG( unsigned )
    void FIX_FLOAT( float )

cdef extern from "platform_dep.h":
    void *alloca(int)

cdef extern from "stdio.h":
    cdef int SEEK_SET
    cdef int SEEK_CUR
    cdef int SEEK_END
    int fseek(FILE *stream, long offset, int whence)
    size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
    long ftell(FILE *stream)
    char *fgets(char *s, int size, FILE *stream)

@cython.boundscheck(False)
@cython.wraparound(False)
def read_and_seek(char *filename, np.int64_t offset1,
                  np.int64_t offset2, np.ndarray buffer, int bytes):
    cdef FILE *f = fopen(filename, "rb")
    cdef void *buf = <void *> buffer.data
    cdef char line[1024]
    cdef size_t n = 1023
    fseek(f, offset1, SEEK_SET)
    fgets(line, n, f)
    fseek(f, offset2, SEEK_CUR)
    fread(buf, 1, bytes, f)
    fclose(f)

def count_art_octs(char *fn, long offset,
                   int min_level, int max_level,
                   int nhydro_vars,
                   level_info):
    cdef int nchild = 8
    cdef int next_record = -1, nLevel = -1
    cdef int dummy_records[9]
    cdef int readin = -1
    cdef FILE *f = fopen(fn, "rb")
    fseek(f, offset, SEEK_SET)
    for _ in range(min_level + 1, max_level + 1):
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
                  np.ndarray[np.int64_t, ndim=2] oct_info):
    #             np.ndarray[np.int64_t, ndim=1] oct_mask,
    #             np.ndarray[np.int64_t, ndim=1] oct_parents,

    # This accepts the filename of the ART header and an integer offset that
    # points to the start of the record *following* the reading of iOctFree and
    # nOct.  For those following along at home, we only need to read:
    #   iOctPr, iOctLv
    cdef int nchild = 8
    cdef int iOct, nLevel, ic1
    cdef np.int64_t next_record = -1
    cdef long long child_record
    cdef int iOctPs[3]
    cdef np.int64_t dummy_records[9]
    cdef int readin = -1
    cdef FILE *f = fopen(fn, "rb")
    fseek(f, offset, SEEK_SET)
    cdef int Level = -1
    cdef int * iNOLL = <int *> alloca(sizeof(int)*(max_level-min_level+1))
    cdef int * iHOLL = <int *> alloca(sizeof(int)*(max_level-min_level+1))
    cdef int iOctMax = 0
    level_offsets = [0]
    for _ in range(min_level + 1, max_level + 1):
        fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
        fread(&Level, sizeof(int), 1, f); FIX_LONG(Level)
        fread(&iNOLL[Level], sizeof(int), 1, f); FIX_LONG(iNOLL[Level])
        fread(&iHOLL[Level], sizeof(int), 1, f); FIX_LONG(iHOLL[Level])
        fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
        iOct = iHOLL[Level] - 1
        nLevel = iNOLL[Level]
        #print "Reading Hierarchy for Level", Lev, Level, nLevel, iOct
        #print ftell(f)
        for ic1 in range(nLevel):
            iOctMax = max(iOctMax, iOct)
            #print readin, iOct, nLevel, sizeof(int)
            next_record = ftell(f)
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            assert readin==52
            next_record += readin + sizeof(int)
            fread(iOctPs, sizeof(int), 3, f)
            FIX_LONG(iOctPs[0]); FIX_LONG(iOctPs[1]); FIX_LONG(iOctPs[2])
            oct_indices[iOct, 0] = iOctPs[0]
            oct_indices[iOct, 1] = iOctPs[1]
            oct_indices[iOct, 2] = iOctPs[2]
            oct_info[iOct, 1] = ic1
            #grid_info[iOct, 2] = iOctPr # we don't seem to need this
            fread(dummy_records, sizeof(int), 6, f) # skip Nb
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            #oct_parents[iOct] = readin - 1
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            oct_levels[iOct] = readin
            fread(&iOct, sizeof(int), 1, f); FIX_LONG(iOct);
            iOct -= 1
            assert next_record > 0
            fseek(f, next_record, SEEK_SET)
            fread(&readin, sizeof(int), 1, f); FIX_LONG(readin)
            assert readin==52

        level_offsets.append(ftell(f))

        #skip over the hydro variables
        #find the length of one child section
        #print 'measuring child record ',
        fread(&next_record, sizeof(int), 1, f);
        #print next_record,
        FIX_LONG(next_record)
        #print next_record
        fseek(f,ftell(f)-sizeof(int),SEEK_SET) #rewind
        #This is a sloppy fix; next_record is 64bit
        #and I don't think FIX_LONG(next_record) is working
        #correctly for 64bits
        if next_record > 4294967296L:
            next_record -= 4294967296L
        assert next_record == 56

        #find the length of all of the children section
        child_record = ftell(f) +  (next_record+2*sizeof(int))*nLevel*nchild
        #print 'Skipping over hydro vars', ftell(f), child_record
        fseek(f, child_record, SEEK_SET)

        # for ic1 in range(nLevel * nchild):
        #     fread(&next_record, sizeof(int), 1, f); FIX_LONG(next_record)
        #     fread(&idc, sizeof(int), 1, f); FIX_LONG(idc); idc -= 1 + (128**3)
        #     fread(&cm, sizeof(int), 1, f); FIX_LONG(cm)
        #     #if cm == 0: oct_mask[idc] = 1
        #     #else: total_masked += 1
        #     assert next_record > 0
        #     fseek(f, next_record - sizeof(int), SEEK_CUR)
    fclose(f)
    return level_offsets

def read_art_root_vars(char *fn, long root_grid_offset,
                    int nhydro_vars, int nx, int ny, int nz,
                    int ix, int iy, int iz, fields, var):

    cdef FILE *f = fopen(fn, "rb")
    cdef int j,l, cell_record_size = nhydro_vars * sizeof(float)
    cdef float temp = -1
    l=0
    fseek(f, root_grid_offset, SEEK_SET)
    # Now we seet out the cell we want
    cdef int my_offset = (((iz * ny) + iy) * nx + ix)
    #print cell_record_size, my_offset, ftell(f)
    fseek(f, cell_record_size * my_offset, SEEK_CUR)
    #(((C)*GridDimension[1]+(B))*GridDimension[0]+A)
    for j in range(nhydro_vars):
        fread(&temp, sizeof(float), 1, f);
        if j in fields:
            FIX_FLOAT(temp)
            var[l]=temp
            l+=1
    fclose(f)

cdef void read_art_vars(FILE *f,
                    int min_level, int max_level, int nhydro_vars,
                    int grid_level,long grid_id,long child_offset,
                    fields,
                    np.ndarray[np.int64_t, ndim=1] level_offsets,
                    var):
    # nhydro_vars is the number of columns- 3 (adjusting for vars)
    # this is normally 10=(8+2chem species)
    cdef int record_size = 2+1+1+nhydro_vars+2
    cdef float temp = -1.0
    cdef float varpad[2]
    cdef int new_padding = -1
    cdef int padding[3]
    cdef long offset = 8*grid_id*record_size*sizeof(float)
    fseek(f, level_offsets[grid_level] + offset, SEEK_SET)
    for j in range(8): #iterate over the children
        l = 0
        fread(padding, sizeof(int), 3, f); FIX_LONG(padding[0])
        #print "Record Size", padding[0]
        # This should be replaced by an fread of nhydro_vars length
        for k in range(nhydro_vars): #iterate over the record
            fread(&temp, sizeof(float), 1, f); FIX_FLOAT(temp)
            #print k, temp
            if k in fields:
                var[j,l] = temp
                l += 1
        fread(varpad, sizeof(float), 2, f)
        fread(&new_padding, sizeof(int), 1, f); FIX_LONG(new_padding)
        assert(padding[0] == new_padding)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def read_art_grid(int varindex,
              np.ndarray[np.int64_t, ndim=1] start_index,
              np.ndarray[np.int32_t, ndim=1] grid_dims,
              np.ndarray[np.float32_t, ndim=3] data,
              np.ndarray[np.uint8_t, ndim=3] filled,
              np.ndarray[np.float32_t, ndim=2] level_data,
              int level, int ref_factor,
              component_grid_info):
    cdef int gi, i, j, k, grid_id
    cdef int ir, jr, kr
    cdef int offi, offj, offk, odind
    cdef np.int64_t di, dj, dk
    cdef np.ndarray[np.int64_t, ndim=1] ogrid_info
    cdef np.ndarray[np.int64_t, ndim=1] og_start_index
    cdef np.float64_t temp_data
    cdef np.int64_t end_index[3]
    cdef int kr_offset, jr_offset, ir_offset
    cdef int to_fill = 0
    # Note that indexing into a cell is:
    #   (k*2 + j)*2 + i
    for i in range(3):
        end_index[i] = start_index[i] + grid_dims[i]
    for gi in range(len(component_grid_info)):
        ogrid_info = component_grid_info[gi]
        grid_id = ogrid_info[1]
        og_start_index = ogrid_info[3:6] #the oct left edge
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
                    if level > 0:
                        odind = (kr*2 + jr)*2 + ir
                        # Replace with an ART-specific reader
                        #temp_data = local_hydro_data.m_var_array[
                        #        level][8*offset + odind]
                        temp_data = level_data[varindex, 8*grid_id + odind]
                    else:
                        kr_offset = kr + <int> (start_index[0] / ref_factor)
                        jr_offset = jr + <int> (start_index[1] / ref_factor)
                        ir_offset = ir + <int> (start_index[2] / ref_factor)
                        odind = (kr_offset * grid_dims[0] + jr_offset)*grid_dims[1] + ir_offset
                        temp_data = level_data[varindex, odind]
                    data[offi, offj, offk] = temp_data
                    filled[offi, offj, offk] = 1
                    to_fill += 1
    return to_fill

@cython.cdivision(True)
@cython.boundscheck(True)
@cython.wraparound(False)
def fill_child_mask(np.ndarray[np.int64_t, ndim=2] file_locations,
                    np.ndarray[np.int64_t, ndim=1] grid_le,
                    np.ndarray[np.uint8_t, ndim=4] art_child_masks,
                    np.ndarray[np.uint8_t, ndim=3] child_mask):

    #loop over file_locations, for each row exracting the index & LE
    #of the oct we will pull pull from art_child_masks
    #then use the art_child_masks info to fill in child_mask
    cdef int i,ioct,x,y,z
    cdef int nocts = file_locations.shape[0]
    cdef int lex,ley,lez
    for i in range(nocts):
        ioct = file_locations[i,1] #from fortran to python indexing?
        lex = file_locations[i,3] - grid_le[0] #the oct left edge x
        ley = file_locations[i,4] - grid_le[1]
        lez = file_locations[i,5] - grid_le[2]
        for x in range(2):
            for y in range(2):
                for z in range(2):
                    child_mask[lex+x,ley+y,lez+z] = art_child_masks[ioct,x,y,z]


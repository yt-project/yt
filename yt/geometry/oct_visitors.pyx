"""
Oct visitor functions




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport numpy
import numpy
from fp_utils cimport *

# Now some visitor functions

cdef void copy_array_f64(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We should always have global_index less than our source.
    # "last" here tells us the dimensionality of the array.
    if selected == 0: return
    cdef int i
    # There are this many records between "octs"
    cdef np.int64_t index = (data.global_index * data.nz)*data.dims
    cdef np.float64_t **p = <np.float64_t**> data.array
    index += oind(data)*data.dims
    for i in range(data.dims):
        p[1][data.index + i] = p[0][index + i]
    data.index += data.dims

cdef void copy_array_i64(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We should always have global_index less than our source.
    # "last" here tells us the dimensionality of the array.
    if selected == 0: return
    cdef int i
    cdef np.int64_t index = (data.global_index * data.nz)*data.dims
    cdef np.int64_t **p = <np.int64_t**> data.array
    index += oind(data)*data.dims
    for i in range(data.dims):
        p[1][data.index + i] = p[0][index + i]
    data.index += data.dims

cdef void count_total_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Count even if not selected.
    # Number of *octs* visited.
    if data.last != o.domain_ind:
        data.index += 1
        data.last = o.domain_ind

cdef void count_total_cells(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Number of *cells* visited and selected.
    data.index += selected

cdef void mark_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We mark them even if they are not selected
    cdef int i
    cdef np.uint8_t *arr = <np.uint8_t *> data.array
    if data.last != o.domain_ind:
        data.last = o.domain_ind
        data.index += 1
    cdef np.int64_t index = data.index * data.nz
    index += oind(data)
    arr[index] = 1

cdef void mask_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    if selected == 0: return
    cdef int i
    cdef np.uint8_t *arr = <np.uint8_t *> data.array
    cdef np.int64_t index = data.global_index * data.nz
    index += oind(data)
    arr[index] = 1

cdef void index_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Note that we provide an index even if the cell is not selected.
    cdef int i
    cdef np.int64_t *arr
    if data.last != o.domain_ind:
        data.last = o.domain_ind
        arr = <np.int64_t *> data.array
        arr[o.domain_ind] = data.index
        data.index += 1

cdef void icoords_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    if selected == 0: return
    cdef np.int64_t *coords = <np.int64_t*> data.array
    cdef int i
    for i in range(3):
        coords[data.index * 3 + i] = (data.pos[i] << data.oref) + data.ind[i]
    data.index += 1

cdef void ires_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    if selected == 0: return
    cdef np.int64_t *ires = <np.int64_t*> data.array
    ires[data.index] = data.level
    data.index += 1

@cython.cdivision(True)
cdef void fcoords_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Note that this does not actually give the correct floating point
    # coordinates.  It gives them in some unit system where the domain is 1.0
    # in all directions, and assumes that they will be scaled later.
    if selected == 0: return
    cdef np.float64_t *fcoords = <np.float64_t*> data.array
    cdef int i
    cdef np.float64_t c, dx 
    dx = 1.0 / ((1 << data.oref) << data.level)
    for i in range(3):
        c = <np.float64_t> ((data.pos[i] << data.oref ) + data.ind[i]) 
        fcoords[data.index * 3 + i] = (c + 0.5) * dx
    data.index += 1

@cython.cdivision(True)
cdef void fwidth_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Note that this does not actually give the correct floating point
    # coordinates.  It gives them in some unit system where the domain is 1.0
    # in all directions, and assumes that they will be scaled later.
    if selected == 0: return
    cdef np.float64_t *fwidth = <np.float64_t*> data.array
    cdef int i
    cdef np.float64_t dx 
    dx = 1.0 / ((1 << data.oref) << data.level)
    for i in range(3):
        fwidth[data.index * 3 + i] = dx
    data.index += 1

cdef void identify_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We assume that our domain has *already* been selected by, which means
    # we'll get all cells within the domain for a by-domain selector and all
    # cells within the domain *and* selector for the selector itself.
    if selected == 0: return
    cdef np.uint8_t *arr = <np.uint8_t *> data.array
    arr[o.domain - 1] = 1

cdef void assign_domain_ind(Oct *o, OctVisitorData *data, np.uint8_t selected):
    o.domain_ind = data.global_index
    data.index += 1

cdef void fill_file_indices_oind(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We fill these arrays, then inside the level filler we use these as
    # indices as we fill a second array from the data.
    if selected == 0: return
    cdef void **p = <void **> data.array
    cdef np.uint8_t *level_arr = <np.uint8_t *> p[0]
    cdef np.int64_t *find_arr = <np.int64_t *> p[1]
    cdef np.uint8_t *cell_arr = <np.uint8_t *> p[2]
    level_arr[data.index] = data.level
    find_arr[data.index] = o.file_ind
    cell_arr[data.index] = oind(data)
    data.index +=1

cdef void fill_file_indices_rind(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We fill these arrays, then inside the level filler we use these as
    # indices as we fill a second array from the data.
    if selected == 0: return
    cdef void **p = <void **> data.array
    cdef np.uint8_t *level_arr = <np.uint8_t *> p[0]
    cdef np.int64_t *find_arr = <np.int64_t *> p[1]
    cdef np.uint8_t *cell_arr = <np.uint8_t *> p[2]
    level_arr[data.index] = data.level
    find_arr[data.index] = o.file_ind
    cell_arr[data.index] = rind(data)
    data.index +=1

cdef void count_by_domain(Oct *o, OctVisitorData *data, np.uint8_t selected):
    cdef np.int64_t *arr
    if selected == 0: return
    # NOTE: We do this for every *cell*.
    arr = <np.int64_t *> data.array
    arr[o.domain - 1] += 1

cdef void store_octree(Oct *o, OctVisitorData *data, np.uint8_t selected):
    cdef np.uint8_t *arr
    if data.last != o.domain_ind:
        data.last = o.domain_ind
        arr = <np.uint8_t *> data.array
        if o.children == NULL:
            arr[data.index] = 0
        if o.children != NULL:
            arr[data.index] = 1
        data.index += 1

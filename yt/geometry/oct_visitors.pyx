"""
Oct visitor functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Author: Christopher Moody <chris.e.moody@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

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

cimport cython
cimport numpy
import numpy
from fp_utils cimport *

# Now some visitor functions

cdef void copy_array_f64(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We should always have global_index less than our source.
    # "last" here tells us the dimensionality of the array.
    if selected == 0: return
    if data.domain > 0 and o.domain != data.domain: return
    cdef int i
    # There are this many records between "octs"
    cdef np.int64_t index = (data.global_index * 8)*data.dims
    cdef np.float64_t **p = <np.float64_t**> data.array
    index += oind(data)*data.dims
    for i in range(data.dims):
        p[1][data.index + i] = p[0][index + i]
    data.index += data.dims

cdef void copy_array_i64(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We should always have global_index less than our source.
    # "last" here tells us the dimensionality of the array.
    if selected == 0: return
    if data.domain > 0 and o.domain != data.domain: return
    cdef int i
    cdef np.int64_t index = (data.global_index * 8)*data.dims
    cdef np.int64_t **p = <np.int64_t**> data.array
    index += oind(data)*data.dims
    for i in range(data.dims):
        p[1][data.index + i] = p[0][index + i]
    data.index += data.dims

cdef void count_total_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Count even if not selected.
    # Number of *octs* visited.
    if data.domain > 0 and o.domain != data.domain: return
    if data.last != o.domain_ind:
        data.index += 1
        data.last = o.domain_ind

cdef void count_total_cells(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Count even if not selected.
    # Number of *octs* visited.
    if data.domain > 0 and o.domain != data.domain: return
    data.index += selected

cdef void mark_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # We mark them even if they are not selected
    cdef int i
    cdef np.uint8_t *arr = <np.uint8_t *> data.array
    if data.last != o.domain_ind:
        data.last = o.domain_ind
        data.index += 1
    cdef np.int64_t index = data.index * 8
    index += oind(data)
    arr[index] = 1

cdef void mask_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    if selected == 0: return
    cdef int i
    cdef np.uint8_t *arr = <np.uint8_t *> data.array
    cdef np.int64_t index = data.global_index * 8
    index += oind(data)
    arr[index] = 1

cdef void index_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    # Note that we provide an index even if the cell is not selected.
    cdef int i
    cdef np.int64_t *arr
    if data.domain > 0 and data.domain != o.domain: return
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
        coords[data.index * 3 + i] = (o.pos[i] << 1) + data.ind[i]
    data.index += 1

cdef void ires_octs(Oct *o, OctVisitorData *data, np.uint8_t selected):
    if selected == 0: return
    cdef np.int64_t *ires = <np.int64_t*> data.array
    ires[data.index] = o.level
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
    dx = 1.0 / (2 << o.level)
    for i in range(3):
        c = <np.float64_t> ((o.pos[i] << 1 ) + data.ind[i]) 
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
    dx = 1.0 / (2 << o.level)
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

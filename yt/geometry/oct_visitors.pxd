"""
Oct visitor definitions file

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
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

cimport numpy as np

cdef struct Oct
cdef struct Oct:
    np.int64_t file_ind     # index with respect to the order in which it was
                            # added
    np.int64_t domain_ind   # index within the global set of domains
    np.int64_t domain       # (opt) addl int index
    Oct **children          # Up to 8 long

cdef struct OctVisitorData:
    np.uint64_t index
    np.uint64_t last
    np.int64_t global_index
    np.int64_t pos[3]       # position in ints
    np.uint8_t ind[3]              # cell position
    void *array
    int dims
    np.int32_t domain
    np.int8_t level

ctypedef void oct_visitor_function(Oct *, OctVisitorData *visitor,
                                   np.uint8_t selected)

cdef oct_visitor_function count_total_octs
cdef oct_visitor_function count_total_cells
cdef oct_visitor_function index_octs
cdef oct_visitor_function icoords_octs
cdef oct_visitor_function ires_octs
cdef oct_visitor_function fcoords_octs
cdef oct_visitor_function fwidth_octs
cdef oct_visitor_function copy_array_f64
cdef oct_visitor_function copy_array_i64
cdef oct_visitor_function identify_octs
cdef oct_visitor_function assign_domain_ind
cdef oct_visitor_function fill_file_indices

cdef inline int cind(int i, int j, int k):
    return (((i*2)+j)*2+k)

cdef inline int oind(OctVisitorData *data):
    return (((data.ind[0]*2)+data.ind[1])*2+data.ind[2])

cdef inline int rind(OctVisitorData *data):
    return (((data.ind[2]*2)+data.ind[1])*2+data.ind[0])

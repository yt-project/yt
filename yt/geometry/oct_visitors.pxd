"""
Oct visitor definitions file

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
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

cimport numpy as np
from selection_routines cimport \
    OctVisitorData, oct_visitor_function
from oct_container cimport \
    Oct

cdef oct_visitor_function count_total_octs
cdef oct_visitor_function count_total_cells
cdef oct_visitor_function mark_octs
cdef oct_visitor_function mask_octs
cdef oct_visitor_function index_octs
cdef oct_visitor_function icoords_octs
cdef oct_visitor_function ires_octs
cdef oct_visitor_function fcoords_octs
cdef oct_visitor_function fwidth_octs
cdef oct_visitor_function copy_array_f64
cdef oct_visitor_function copy_array_i64

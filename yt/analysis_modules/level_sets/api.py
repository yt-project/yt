"""
API for level_sets

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: J.S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
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

from .contour_finder import \
    GridConsiderationQueue, \
    coalesce_join_tree, \
    identify_contours

from .clump_handling import \
    Clump, \
    find_clumps, \
    get_lowest_clumps, \
    write_clump_hierarchy, \
    write_clumps, \
    write_old_clump_hierarchy, \
    write_old_clumps, \
    write_old_clump_info, \
    _DistanceToMainClump

from .clump_tools import \
    recursive_all_clumps, \
    return_all_clumps, \
    return_bottom_clumps, \
    recursive_bottom_clumps, \
    clump_list_sort

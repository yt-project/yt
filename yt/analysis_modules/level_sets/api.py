"""
API for level_sets



"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import issue_deprecation_warning

issue_deprecation_warning(
    "The level_sets module has been moved to yt.data_objects.level_sets."
    "This import is deprecated and will be removed in a future release."
    "Please, change the import in your scripts from "
    "'from yt.analysis_modules.level_sets' to "
    "'from yt.data_objects.level_sets.'.")

from yt.data_objects.level_sets.contour_finder import \
    identify_contours

from yt.data_objects.level_sets.clump_handling import \
    Clump, \
    find_clumps, \
    get_lowest_clumps

from yt.data_objects.level_sets.clump_info_items import \
    add_clump_info

from yt.data_objects.level_sets.clump_validators import \
    add_validator
    
from yt.data_objects.level_sets.clump_tools import \
    recursive_all_clumps, \
    return_all_clumps, \
    return_bottom_clumps, \
    recursive_bottom_clumps, \
    clump_list_sort

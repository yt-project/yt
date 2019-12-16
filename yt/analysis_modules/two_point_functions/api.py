"""
API for two_point_functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2020, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


from yt.funcs import issue_deprecation_warning

issue_deprecation_warning(
    "The two_point_functions module has been deprecated. This code has been "
    "moved to the yt attic (https://github.com/yt-project/yt_attic) and "
    "will be removed in a future release.")

from .two_point_functions import \
    TwoPointFunctions, \
    FcnSet

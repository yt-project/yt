"""
API for halo_finding



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import issue_deprecation_warning

issue_deprecation_warning(
    "Development of the halo_finding module has been moved to "
    "the yt_astro_analysis package. This version is deprecated "
    "and will be removed from yt in a future release. See "
    "https://github.com/yt-project/yt_astro_analysis for further "
    "information.")

from .halo_objects import \
    Halo, \
    HOPHalo, \
    LoadedHalo, \
    FOFHalo, \
    HaloList, \
    HOPHaloList, \
    FOFHaloList, \
    LoadedHaloList, \
    GenericHaloFinder, \
    HOPHaloFinder, \
    FOFHaloFinder, \
    HaloFinder, \
    LoadHaloes, \
    LoadTextHalos, \
    LoadTextHaloes

"""
API for halo_analysis



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import issue_deprecation_warning

issue_deprecation_warning(
    "Development of the HaloCatalog module has been moved to "
    "the yt_astro_analysis package. This version is deprecated "
    "and will be removed from yt in a future release. See "
    "https://github.com/yt-project/yt_astro_analysis for further "
    "information.")

from .halo_catalog import \
    HaloCatalog

from .halo_callbacks import \
    add_callback

from .halo_finding_methods import \
    add_finding_method

from .halo_filters import \
    add_filter
     
from .halo_quantities import \
    add_quantity

from .halo_recipes import \
    add_recipe

"""
API for cosmology analysis.



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
    "Development of the CosmologySplice module has been moved to "
    "the yt_astro_analysis package. This version is deprecated "
    "and will be removed from yt in a future release. See "
    "https://github.com/yt-project/yt_astro_analysis for further "
    "information.")

from .cosmology_splice import \
    CosmologySplice

from .light_cone.api import \
    LightCone

from .light_ray.api import \
    LightRay

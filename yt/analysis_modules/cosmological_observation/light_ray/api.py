"""
API for light_ray



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
    "Development of the LightRay module has been moved to the Trident "
    "package. This version is deprecated and will be removed from yt "
    "in a future release. See https://github.com/trident-project/trident "
    "for further information.")

from .light_ray import \
    LightRay

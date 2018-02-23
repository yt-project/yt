"""
API for halo_mass_function



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
    "The halo_mass_function module does not function correctly and has been "
    "deprecated. This code has been moved to the yt attic "
    "(https://github.com/yt-project/yt_attic) and will be removed in a "
    "future release.")

from .halo_mass_function import \
    HaloMassFcn, \
    TransferFunction, \
    integrate_inf

"""
API for yt.frontends.fits
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .data_structures import \
      FITSGrid, \
      FITSHierarchy, \
      FITSStaticOutput

from .fields import \
      FITSFieldInfo, \
      add_fits_field

from .io import \
      IOHandlerFITS

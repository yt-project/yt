"""
API for yt.frontends.gdf



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from .data_structures import \
      GDFGrid, \
      GDFHierarchy, \
      GDFDataset

from .fields import \
      GDFFieldInfo
add_gdf_field = GDFFieldInfo.add_field

from .io import \
      IOHandlerGDFHDF5

from . import tests

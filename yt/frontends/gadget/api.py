"""
API for yt.frontends.gadget



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .data_structures import \
      GadgetGrid, \
      GadgetHierarchy, \
      GadgetStaticOutput

from .fields import \
      GadgetFieldInfo, \
      add_gadget_field

from .io import \
      IOHandlerGadget

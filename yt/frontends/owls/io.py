"""
OWLS data-file handling function




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.gadget.io import \
    IOHandlerGadgetHDF5

class IOHandlerOWLS(IOHandlerGadgetHDF5):
    _dataset_type = "OWLS"

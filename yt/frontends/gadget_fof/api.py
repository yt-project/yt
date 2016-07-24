"""
API for GadgetFOF frontend




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .data_structures import \
    GadgetFOFParticleIndex, \
    GadgetFOFHDF5File, \
    GadgetFOFDataset, \
    GadgetFOFHaloParticleIndex, \
    GadgetFOFHaloDataset, \
    GagdetFOFHaloContainer

from .io import \
    IOHandlerGadgetFOFHDF5, \
    IOHandlerGadgetFOFHaloHDF5

from .fields import \
    GadgetFOFFieldInfo, \
    GadgetFOFHaloFieldInfo

from . import tests

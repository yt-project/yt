"""
API for gadget_fof_plus frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Britton Smith.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import \
    mylog

mylog.warn("GizmoDataset dataset overrides GadgetHDF5Dataset.")

from .data_structures import \
    GizmoDataset

from .fields import \
    GizmoFieldInfo

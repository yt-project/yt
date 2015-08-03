"""
Some old field names.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.units.yt_array import YTArray
from yt.units.unit_object import Unit
from yt.utilities.exceptions import YTUnitOperationError


class Quantity(YTArray):
    """
    A physical quantity. Attaches units to a scalar.

    """
    def __new__(cls, input_array, input_units=None):
        if isinstance(input_array, Quantity):
            return input_array

        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)

        # Restrict the array to a scalar.
        if obj.size != 1:
            raise ValueError("A Quantity can contain only one element. The "
                "caller provided the array %s with %s elements."
                % (obj, obj.size))

        return YTArray.__new__(cls, input_array, input_units)

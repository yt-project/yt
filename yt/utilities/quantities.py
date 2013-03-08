"""
Quantities -- floats with units.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley

Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Casey W. Stark.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import numpy as np

from yt.data_objects.yt_array import YTArray
from yt.utilities.units import Unit, UnitOperationError


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

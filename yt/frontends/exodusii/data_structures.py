"""
Data structures for Exodusii.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import string
from itertools import takewhile
from collections import OrderedDict
import numpy as np
from yt.funcs import *
from netCDF4 import Dataset
from yt.data_objects.static_output import \
           Dataset
from .fields import ExodusiiFieldInfo


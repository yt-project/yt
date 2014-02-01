"""
General field-related functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

def get_radius(data, field_prefix):
    center = data.get_field_parameter("center").in_units("cm")
    DW = (data.pf.domain_right_edge - data.pf.domain_left_edge).in_units("cm")
    radius = data.pf.arr(np.zeros(data[field_prefix+"x"].shape,
                         dtype='float64'), 'cm')
    r = radius.copy()
    if any(data.pf.periodicity):
        rdw = radius.copy()
    for i, ax in enumerate('xyz'):
        np.subtract(data["%s%s" % (field_prefix, ax)].in_units("cm"),
                    center[i], r)
        if data.pf.periodicity[i] == True:
            np.abs(r, r)
            np.subtract(r, DW[i], rdw)
            np.abs(rdw, rdw)
            np.minimum(r, rdw, r)
        np.power(r, 2.0, r)
        np.add(radius, r, radius)
        if data.pf.dimensionality < i+1:
            break
    np.sqrt(radius, radius)
    return radius

"""
Fields from interpolating data tables.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.local_fields import add_field

from yt.utilities.linear_interpolators import \
    UnilinearFieldInterpolator, \
    BilinearFieldInterpolator, \
    TrilinearFieldInterpolator

_int_class = {1: UnilinearFieldInterpolator,
              2: BilinearFieldInterpolator,
              3: TrilinearFieldInterpolator}

def add_interpolated_field(name, units, table_data, axes_data, axes_fields,
                           ftype="gas", particle_type=False, validators=None,
                           truncate=True):
    
    if len(table_data.shape) not in _int_class:
        raise RuntimeError("Interpolated field can only be created from 1d, 2d, or 3d data.")

    if len(axes_fields) != len(axes_data) or len(axes_fields) != len(table_data.shape):
        raise RuntimeError("Data dimension mismatch: data is %d, %d axes data provided, and %d axes fields provided." %
                           (len(table_data.shape), len(axes_data), len(axes_fields)))

    int_class = _int_class[len(table_data.shape)]
    my_interpolator = int_class(table_data, axes_data, axes_fields, truncate=truncate)

    def _interpolated_field(field, data):
        return my_interpolator(data)
    add_field((ftype, name), 
              function=_interpolated_field,
              units=units,
              validators=validators, particle_type=particle_type)
              

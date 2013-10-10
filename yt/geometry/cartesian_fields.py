"""
Cartesian fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    ValidateDataField, \
    ValidateGridType, \
    ValidateParameter, \
    ValidateSpatial, \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter
from yt.data_objects.yt_array import YTArray

CartesianFieldInfo = FieldInfoContainer()
CartesianFieldInfo.name = id(CartesianFieldInfo)
add_cart_field = CartesianFieldInfo.add_field

def _dx(field, data):
    return YTArray(data.fwidth[...,0], 'code_length',
                   registry = data.pf.unit_registry)
add_cart_field('dx', function=_dx, display_field=False, units='code_length')

def _dy(field, data):
    return YTArray(data.fwidth[...,1], 'code_length',
                   registry = data.pf.unit_registry)
add_cart_field('dy', function=_dy, display_field=False, units='code_length')

def _dz(field, data):
    return YTArray(data.fwidth[...,2], 'code_length',
                   registry = data.pf.unit_registry)
add_cart_field('dz', function=_dz, display_field=False, units='code_length')

def _coordX(field, data):
    return YTArray(data.fcoords[...,0], 'code_length',
                   registry = data.pf.unit_registry)
add_cart_field('x', function=_coordX, display_field=False, units='code_length')

def _coordY(field, data):
    return YTArray(data.fcoords[...,1], 'code_length',
                   registry = data.pf.unit_registry)
add_cart_field('y', function=_coordY, display_field=False, units='code_length')

def _coordZ(field, data):
    return YTArray(data.fcoords[...,2], 'code_length',
                   registry = data.pf.unit_registry)
add_cart_field('z', function=_coordZ, display_field=False, units='code_length')


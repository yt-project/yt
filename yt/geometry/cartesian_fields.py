"""
Cartesian fields

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

CartesianFieldInfo = FieldInfoContainer()
CartesianFieldInfo.name = id(CartesianFieldInfo)
add_cart_field = CartesianFieldInfo.add_field

def _dx(field, data):
    return data.pf.domain_width[0] * data.fwidth[...,0]
add_cart_field('dx', function=_dx, display_field=False)

def _dy(field, data):
    return data.pf.domain_width[1] * data.fwidth[...,1]
add_cart_field('dy', function=_dy, display_field=False)

def _dz(field, data):
    return data.pf.domain_width[2] * data.fwidth[...,2]
add_cart_field('dz', function=_dz, display_field=False)

def _coordX(field, data):
    return data.pf.domain_left_edge[0] + data.fcoords[...,0]
add_cart_field('x', function=_coordX, display_field=False)

def _coordY(field, data):
    return data.pf.domain_left_edge[1] + data.fcoords[...,1]
add_cart_field('y', function=_coordY, display_field=False)

def _coordZ(field, data):
    return data.pf.domain_left_edge[2] + data.fcoords[...,2]
add_cart_field('z', function=_coordZ, display_field=False)


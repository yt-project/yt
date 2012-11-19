"""
Cylindrical fields

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

CylindricalFieldInfo = FieldInfoContainer()
CylindricalFieldInfo.name = id(CylindricalFieldInfo)
add_cyl_field = CylindricalFieldInfo.add_field

#
# Cylindrical fields
#

def _unknown_coord(field, data):
    raise YTCoordinateNotImplemented
add_cyl_field("dx", function=_unknown_coord)
add_cyl_field("dy", function=_unknown_coord)
add_cyl_field("x", function=_unknown_coord)
add_cyl_field("y", function=_unknown_coord)

def _dr(field, data):
    return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[0]
add_cyl_field('dr', function=_dr, display_field=False,
          validators=[ValidateSpatial(0)])

def _dz(field, data):
    return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[1]
add_cyl_field('dz', function=_dz,
          display_field=False, validators=[ValidateSpatial(0)])

def _dtheta(field, data):
    return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[2]
add_cyl_field('dtheta', function=_dtheta,
          display_field=False, validators=[ValidateSpatial(0)])

def _coordR(field, data):
    dim = data.ActiveDimensions[0]
    return (np.ones(data.ActiveDimensions, dtype='float64')
                   * np.arange(data.ActiveDimensions[0])[:,None,None]
            +0.5) * data['dr'] + data.LeftEdge[0]
add_cyl_field('r', function=_coordR, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordZ(field, data):
    dim = data.ActiveDimensions[1]
    return (np.ones(data.ActiveDimensions, dtype='float64')
                   * np.arange(data.ActiveDimensions[1])[None,:,None]
            +0.5) * data['dz'] + data.LeftEdge[1]
add_cyl_field('z', function=_coordZ, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordTheta(field, data):
    dim = data.ActiveDimensions[2]
    return (np.ones(data.ActiveDimensions, dtype='float64')
                   * np.arange(data.ActiveDimensions[2])[None,None,:]
            +0.5) * data['dtheta'] + data.LeftEdge[2]
add_cyl_field('theta', function=_coordTheta, display_field=False,
          validators=[ValidateSpatial(0)])

def _CylindricalVolume(field, data):
    return data["dtheta"] * data["r"] * data["dr"] * data["dz"]
add_cyl_field("CellVolume", function=_CylindricalVolume)


## Now the Polar fields

PolarFieldInfo = FieldInfoContainer()
PolarFieldInfo.name = id(PolarFieldInfo)
add_pol_field = PolarFieldInfo.add_field

add_pol_field("dx", function=_unknown_coord)
add_pol_field("dy", function=_unknown_coord)
add_pol_field("x", function=_unknown_coord)
add_pol_field("y", function=_unknown_coord)

def _dr(field, data):
    return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[0]
add_pol_field('dr', function=_dr, display_field=False,
          validators=[ValidateSpatial(0)])

def _dtheta(field, data):
    return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[1]
add_pol_field('dtheta', function=_dtheta,
          display_field=False, validators=[ValidateSpatial(0)])

def _coordR(field, data):
    dim = data.ActiveDimensions[0]
    return (np.ones(data.ActiveDimensions, dtype='float64')
                   * np.arange(data.ActiveDimensions[0])[:,None,None]
            +0.5) * data['dr'] + data.LeftEdge[0]
add_pol_field('r', function=_coordR, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordTheta(field, data):
    dim = data.ActiveDimensions[2]
    return (np.ones(data.ActiveDimensions, dtype='float64')
                   * np.arange(data.ActiveDimensions[1])[None,:,None]
            +0.5) * data['dtheta'] + data.LeftEdge[1]
add_pol_field('theta', function=_coordTheta, display_field=False,
          validators=[ValidateSpatial(0)])

def _CylindricalVolume(field, data):
    return data["dtheta"] * data["r"] * data["dr"] * data["dz"]
add_pol_field("CellVolume", function=_CylindricalVolume)

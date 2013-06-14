"""
Pluto-specific fields

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2009-2011 J. S. Oishi, Matthew Turk.  All Rights Reserved.

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

from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    FieldInfo, \
    NullFunc, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
import numpy as np

KnownPlutoFields = FieldInfoContainer()
add_pluto_field = KnownPlutoFields.add_field

PlutoFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = PlutoFieldInfo.add_field

add_pluto_field("rho", function=NullFunc, take_log=True,
                 validators = [ValidateDataField("density")],
                 units=r"\rm{g}/\rm{cm}^3")

KnownPlutoFields["rho"]._projected_units =r"\rm{g}/\rm{cm}^2"

add_pluto_field("vx1", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("X-Momentum")],
                 units=r"",display_name=r"M_x")
KnownPlutoFields["vx1"]._projected_units=r""

add_pluto_field("vx2", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("Y-Momentum")],
                 units=r"",display_name=r"M_y")
KnownPlutoFields["vx2"]._projected_units=r""

add_pluto_field("vx3", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("Z-Momentum")],
                 units=r"",display_name=r"M_z")
KnownPlutoFields["vx3"]._projected_units=r""

add_pluto_field("prs", function=NullFunc, take_log=True,
                 validators = [ValidateDataField("energy-density")],
                 units=r"\rm{erg}/\rm{cm}^3")

def _Density(field,data):
    """A duplicate of the density field. This is needed because when you try 
    to instantiate a PlotCollection without passing in a center, the code
    will try to generate one for you using the "Density" field, which gives an error 
    if it isn't defined.

    """
    return data["rho"]
add_field("Density",function=_Density, take_log=True,
          units=r'\rm{g}/\rm{cm^3}')

def _Xmomentum(field, data):
    """ Generate x-momentum. """
    return data["vx1"]*data["density"]
add_field("X-momentum",function=_Xmomentum, take_log=False,
          units=r'\rm{g}/\rm{cm^2 s}')

def _Ymomentum(field, data):
    """ Generate y-momentum  """
    return data["vx2"]*data["density"]
add_field("Y-momentum",function=_Ymomentum, take_log=False,
          units=r'\rm{g}/\rm{cm^2 s}')

def _Zmomentum(field,data):
    """ Generate z-momentum"""
    return data["vx3"]*data["density"]
add_field("Z-Momentum",function=_Zmomentum, take_log=False,
          units=r'\rm{g}/\rm{cm^2 s}')


"""
Chombo-specific fields

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

KnownChomboFields = FieldInfoContainer()
add_chombo_field = KnownChomboFields.add_field

ChomboFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = ChomboFieldInfo.add_field

add_chombo_field("density", function=NullFunc, take_log=True,
                 validators = [ValidateDataField("density")],
                 units="g/cm**3")

add_chombo_field("X-momentum", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("X-Momentum")],
                 units="g/cm**2/s",display_name=r"M_x")

add_chombo_field("Y-momentum", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("Y-Momentum")],
                 units="g/cm**2/s",display_name=r"M_y")

add_chombo_field("Z-momentum", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("Z-Momentum")],
                 units="g/cm**2/s",display_name=r"M_z")

add_chombo_field("X-magnfield", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("X-Magnfield")],
                 units="gauss",display_name=r"B_x")

add_chombo_field("Y-magnfield", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("Y-Magnfield")],
                 units="gauss",display_name=r"B_y")

add_chombo_field("Z-magnfield", function=NullFunc, take_log=False,
                  validators = [ValidateDataField("Z-Magnfield")],
                  units="gauss",display_name=r"B_z")

add_chombo_field("energy-density", function=NullFunc, take_log=True,
                 validators = [ValidateDataField("energy-density")],
                 units="erg/cm**3")

add_chombo_field("radiation-energy-density", function=NullFunc, take_log=True,
                 validators = [ValidateDataField("radiation-energy-density")],
                 units="erg/cm**3")

def _Density(field,data):
    """A duplicate of the density field. This is needed because when you try 
    to instantiate a PlotCollection without passing in a center, the code
    will try to generate one for you using the "Density" field, which gives an error 
    if it isn't defined.

    """
    return data["density"]
add_field("Density",function=_Density, take_log=True,
          units='g/cm**3')

def _Bx(field,data):
    return data["X-magnfield"]
add_field("Bx", function=_Bx, take_log=False,
          units="gauss", display_name=r"B_x")

def _By(field,data):
    return data["Y-magnfield"]
add_field("By", function=_By, take_log=False,
          units="gauss", display_name=r"B_y")

def _Bz(field,data):
    return data["Z-magnfield"]
add_field("Bz", function=_Bz, take_log=False,
          units="gauss", display_name=r"B_z")

def _MagneticEnergy(field,data):
    return (data["X-magnfield"]**2 +
            data["Y-magnfield"]**2 +
            data["Z-magnfield"]**2)/2.
add_field("MagneticEnergy", function=_MagneticEnergy, take_log=True,
          units=r"erg/cm**3", display_name=r"B^2 / 8 \pi")

def _xVelocity(field, data):
    """ Generate x-velocity from x-momentum and density. """
    return data["X-momentum"]/data["density"]
add_field("x-velocity",function=_xVelocity, take_log=False,
          units='cm/s')

def _yVelocity(field,data):
    """ Generate y-velocity from y-momentum and density. """
    #try:
    #    return data["xvel"]
    #except KeyError:
    return data["Y-momentum"]/data["density"]
add_field("y-velocity",function=_yVelocity, take_log=False,
          units='cm/s')

def _zVelocity(field,data):
    """ Generate z-velocity from z-momentum and density. """
    return data["Z-momentum"]/data["density"]
add_field("z-velocity",function=_zVelocity, take_log=False,
          units='cm/s')

def particle_func(p_field, dtype='float64'):
    def _Particles(field, data):
        io = data.hierarchy.io
        if not data.NumberOfParticles > 0:
            return np.array([], dtype=dtype)
        else:
            return io._read_particles(data, p_field).astype(dtype)
        
    return _Particles

_particle_field_list = ["mass",
                        "position_x",
                        "position_y",
                        "position_z",
                        "momentum_x",
                        "momentum_y",
                        "momentum_z",
                        "angmomen_x",
                        "angmomen_y",
                        "angmomen_z",
                        "mlast",
                        "mdeut",
                        "n",
                        "mdot",
                        "burnstate",
                        "id"]

for pf in _particle_field_list:
    pfunc = particle_func("particle_%s" % (pf))
    add_field("particle_%s" % pf, function=pfunc,
              validators = [ValidateSpatial(0)],
              particle_type=True)

"""
Orion-specific fields

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: UC Berkeley
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 J. S. Oishi, Matthew Turk.  All Rights Reserved.

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


from UniversalFields import *

add_field = add_orion_field

add_field("density", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("density")],
          units=r"\rm{g}/\rm{cm}^3")

add_field("eden", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("eden")],
          units=r"\rm{erg}/\rm{cm}^3")

translation_dict = {"x-velocity": "xvel",
                    "y-velocity": "yvel",
                    "z-velocity": "zvel",
                    "Density": "density",
                    "Total_Energy": "eden",
                    "Temperature": "temperature",
                    "x-momentum": "xmom",
                    "y-momentum": "ymom",
                    "z-momentum": "zmom"
                   }

def _generate_translation(mine, theirs):
    add_field(theirs, function=lambda a, b: b[mine], take_log=True)

for f,v in translation_dict.items():
    if v not in OrionFieldInfo:
        add_field(v, function=lambda a,b: None, take_log=True,
                  validators = [ValidateDataField(v)])
    #print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)

def _xVelocity(field, data):
    """generate x-velocity from x-momentum and density

    """
    return data["xmom"]/data["density"]
add_field("x-velocity",function=_xVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

def _yVelocity(field,data):
    """generate y-velocity from y-momentum and density

    """
    #try:
    #    return data["xvel"]
    #except KeyError:
    return data["ymom"]/data["density"]
add_field("y-velocity",function=_yVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

def _zVelocity(field,data):
    """generate z-velocity from z-momentum and density

    """
    return data["zmom"]/data["density"]
add_field("z-velocity",function=_zVelocity, take_log=False,
          units=r'\rm{cm}/\rm{s}')

def _ThermalEnergy(field, data):
    """generate thermal (gas energy). Dual Energy Formalism was
        implemented by Stella, but this isn't how it's called, so I'll
        leave that commented out for now.
    """
    #if data.pf["DualEnergyFormalism"]:
    #    return data["Gas_Energy"]
    #else:
    return data["Total_Energy"] - 0.5 * data["density"] * (
        data["x-velocity"]**2.0
        + data["y-velocity"]**2.0
        + data["z-velocity"]**2.0 )
add_field("ThermalEnergy", function=_ThermalEnergy,
          units=r"\rm{ergs}/\rm{cm^3}")

def _Pressure(field,data):
    """M{(Gamma-1.0)*e, where e is thermal energy density
       NB: this will need to be modified for radiation
    """
    return (data.pf["Gamma"] - 1.0)*data["ThermalEnergy"]
add_field("Pressure", function=_Pressure, units=r"\rm{dyne}/\rm{cm}^{2}")

def _Temperature(field,data):
    return (data.pf["Gamma"]-1.0)*data.pf["mu"]*mh*data["ThermalEnergy"]/(kboltz*data["Density"])
add_field("Temperature",function=_Temperature,units=r"\rm{Kelvin}",take_log=False)

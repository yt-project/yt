"""
Chombo-specific fields

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: UC Berkeley
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 J. S. Oishi, Matthew Turk.  All Rights Reserved.

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
class FLASHFieldContainer(CodeFieldInfoContainer):
    _shared_state = {}
    _field_list = {}
FLASHFieldInfo = FLASHFieldContainer()
add_flash_field = FLASHFieldInfo.add_field

add_field = add_flash_field

# Common fields in FLASH: (Thanks to John ZuHone for this list)
#
# dens gas mass density (g/cc) --
# eint internal energy (ergs/g) --
# ener total energy (ergs/g), with 0.5*v^2 --
# gamc gamma defined as ratio of specific heats, no units
# game gamma defined as in , no units
# gpol gravitational potential from the last timestep (ergs/g)
# gpot gravitational potential from the current timestep (ergs/g)
# grac gravitational acceleration from the current timestep (cm s^-2)
# pden particle mass density (usually dark matter) (g/cc)
# pres pressure (erg/cc)
# temp temperature (K) --
# velx velocity x (cm/s) --
# vely velocity y (cm/s) --
# velz velocity z (cm/s) --

translation_dict = {"x-velocity": "velx",
                    "y-velocity": "vely",
                    "z-velocity": "velz",
                    "Density": "dens",
                    "Total_Energy": "ener",
                    "Gas_Energy": "eint",
                    "Temperature": "temp",
                   }

def _generate_translation(mine, theirs):
    add_field(theirs, function=lambda a, b: b[mine], take_log=True)

for f,v in translation_dict.items():
    if v not in FLASHFieldInfo:
        add_field(v, function=lambda a,b: None, take_log=False,
                  validators = [ValidateDataField(v)])
    #print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)

add_field("gamc", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("gamc")],
          units = r"\rm{ratio\/of\/specific\/heats}")

add_field("game", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("game")],
          units = r"\rm{ratio\/of\/specific\/heats}")

add_field("gpot", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("gpot")],
          units = r"\rm{ergs\//\/g}")

add_field("gpot", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("gpol")],
          units = r"\rm{ergs\//\/g}")

add_field("grac", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("grac")],
          units = r"\rm{cm\/s^{-2}}")

add_field("pden", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("pden")],
          units = r"\rm{g}\//\/\rm{cm}^{3}")

add_field("pres", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("pres")],
          units = r"\rm{erg}\//\/\rm{cm}^{3}")

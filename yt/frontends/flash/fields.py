"""
FLASH-specific fields

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk, John ZuHone.  All Rights Reserved.

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
    CodeFieldInfoContainer, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields

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
                    "particle_position_x" : "particle_posx",
                    "particle_position_y" : "particle_posy",
                    "particle_position_z" : "particle_posz",
                    "particle_velocity_x" : "particle_velx",
                    "particle_velocity_y" : "particle_vely",
                    "particle_velocity_z" : "particle_velz",
                    "particle_index" : "particle_tag",
                    "Electron_Fraction" : "elec",
                    "HI_Fraction" : "h   ",
                    "HD_Fraction" : "hd  ",
                    "HeI_Fraction": "hel ",
                    "HeII_Fraction": "hep ",
                    "HeIII_Fraction": "hepp",
                    "HM_Fraction": "hmin",
                    "HII_Fraction": "hp  ",
                    "H2I_Fraction": "htwo",
                    "H2II_Fraction": "htwp",
                    "DI_Fraction": "deut",
                    "DII_Fraction": "dplu",
                    "ParticleMass": "particle_mass"}

def _get_density(fname):
    def _dens(field, data):
        return data[fname] * data['Density']
    return _dens

for fn1, fn2 in translation_dict.items():
    if fn1.endswith("_Fraction"):
        add_field(fn1.split("_")[0] + "_Density",
                  function=_get_density(fn1), take_log=True)

def _get_alias(alias):
    def _alias(field, data):
        return data[alias]
    return _alias

def _generate_translation(mine, theirs):
    pfield = theirs.startswith("particle")
    add_field(theirs, function=_get_alias(mine), take_log=True,
              particle_type = pfield)

def _get_convert(fname):
    def _conv(data):
        return data.convert(fname)
    return _conv

add_field("dens", function=lambda a,b: None, take_log=True,
          convert_function=_get_convert("dens"),
          units=r"\rm{g}/\rm{cm}^3")
add_field("xvel", function=lambda a,b: None, take_log=False,
          convert_function=_get_convert("xvel"),
          units=r"\rm{cm}/\rm{s}")
add_field("yvel", function=lambda a,b: None, take_log=False,
          convert_function=_get_convert("yvel"),
          units=r"\rm{cm}/\rm{s}")
add_field("zvel", function=lambda a,b: None, take_log=False,
          convert_function=_get_convert("zvel"),
          units=r"\rm{cm}/\rm{s}")
add_field("particle_xvel", function=lambda a,b: None, take_log=False,
          convert_function=_get_convert("particle_xvel"),
          units=r"\rm{cm}/\rm{s}")
add_field("particle_yvel", function=lambda a,b: None, take_log=False,
          convert_function=_get_convert("particle_yvel"),
          units=r"\rm{cm}/\rm{s}")
add_field("particle_zvel", function=lambda a,b: None, take_log=False,
          convert_function=_get_convert("particle_zvel"),
          units=r"\rm{cm}/\rm{s}")
add_field("temp", function=lambda a,b: None, take_log=True,
          convert_function=_get_convert("temp"),
          units=r"\rm{K}")
add_field("pres", function=lambda a,b: None, take_log=True,
          convert_function=_get_convert("pres"),
          units=r"\rm{unknown}")

for f,v in translation_dict.items():
    if v not in FLASHFieldInfo:
        pfield = v.startswith("particle")
        add_field(v, function=lambda a,b: None, take_log=False,
                  validators = [ValidateDataField(v)],
                  particle_type = pfield)
    #print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)

add_field("gamc", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("gamc")],
          units = r"\rm{ratio\/of\/specific\/heats}")

add_field("game", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("game")],
          units = r"\rm{ratio\/of\/specific\/heats}")

add_field("gpot", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("gpot")],
          units = r"\rm{ergs\//\/g}")

add_field("gpol", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("gpol")],
          units = r"\rm{ergs\//\/g}")

add_field("grac", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("grac")],
          units = r"\rm{cm\/s^{-2}}")

add_field("pden", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("pden")],
          units = r"\rm{g}\//\/\rm{cm}^{3}")

add_field("pres", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("pres")],
          units = r"\rm{erg}\//\/\rm{cm}^{3}")

add_field("magx", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("magx")],
          units = r"\rm{G}")

add_field("magy", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("magy")],
          units = r"\rm{G}")

add_field("magz", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("magz")],
          units = r"\rm{G}")

add_field("magp", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("magp")],
          units = r"\rm{erg}\//\/\rm{cm}^{3}")

add_field("divb", function=lambda a,b: None, take_log=False,
          validators = [ValidateDataField("divb")],
          units = r"\rm{G}\/\rm{cm}")

def _convertParticleMassMsun(data):
    return 1.0/1.989e33
def _ParticleMassMsun(field, data):
    return data["ParticleMass"]
add_field("ParticleMassMsun",
          function=_ParticleMassMsun, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMassMsun,
          particle_convert_function=_ParticleMassMsun)

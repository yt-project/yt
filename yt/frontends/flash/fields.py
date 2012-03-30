"""
FLASH-specific fields

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk, John ZuHone.  All Rights Reserved.

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
    NullFunc, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
from yt.utilities.physical_constants import \
    kboltz
KnownFLASHFields = FieldInfoContainer()
add_flash_field = KnownFLASHFields.add_field

FLASHFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = FLASHFieldInfo.add_field

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
                    "Temperature": "temp",
                    "Pressure" : "pres", 
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
                    "ParticleMass": "particle_mass",
                    "Flame_Fraction": "flam"}

def _get_density(fname):
    def _dens(field, data):
        return data[fname] * data['Density']
    return _dens

for fn1, fn2 in translation_dict.items():
    if fn1.endswith("_Fraction"):
        add_field(fn1.split("_")[0] + "_Density",
                  function=_get_density(fn1), take_log=True,
                  display_name="%s\/Density" % fn1.split("_")[0])

def _get_convert(fname):
    def _conv(data):
        return data.convert(fname)
    return _conv

add_flash_field("dens", function=NullFunc, take_log=True,
                convert_function=_get_convert("dens"),
                units=r"\rm{g}/\rm{cm}^3")
add_flash_field("velx", function=NullFunc, take_log=False,
                convert_function=_get_convert("velx"),
                units=r"\rm{cm}/\rm{s}")
add_flash_field("vely", function=NullFunc, take_log=False,
                convert_function=_get_convert("vely"),
                units=r"\rm{cm}/\rm{s}")
add_flash_field("velz", function=NullFunc, take_log=False,
                convert_function=_get_convert("velz"),
                units=r"\rm{cm}/\rm{s}")
add_flash_field("ener", function=NullFunc, take_log=True,
                convert_function=_get_convert("ener"),
                units=r"\rm{erg}/\rm{g}")
add_flash_field("eint", function=NullFunc, take_log=True,
                convert_function=_get_convert("eint"),
                units=r"\rm{erg}/\rm{g}")
add_flash_field("particle_posx", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_posx"),
                units=r"\rm{cm}", particle_type=True)
add_flash_field("particle_posy", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_posy"),
                units=r"\rm{cm}", particle_type=True)
add_flash_field("particle_posz", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_posz"),
                units=r"\rm{cm}", particle_type=True)
add_flash_field("particle_velx", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_velx"),
                units=r"\rm{cm}/\rm{s}", particle_type=True)
add_flash_field("particle_vely", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_vely"),
                units=r"\rm{cm}/\rm{s}", particle_type=True)
add_flash_field("particle_velz", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_velz"),
                units=r"\rm{cm}/\rm{s}", particle_type=True)
add_flash_field("particle_tag", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_tag"),
                particle_type=True)
add_flash_field("particle_mass", function=NullFunc, take_log=False,
                convert_function=_get_convert("particle_mass"),
                units=r"\rm{g}", particle_type=True)
add_flash_field("temp", function=NullFunc, take_log=True,
                convert_function=_get_convert("temp"),
                units=r"\rm{K}")
add_flash_field("tele", function=NullFunc, take_log=True,
                convert_function=_get_convert("tele"),
                units = r"\rm{K}")
add_flash_field("pres", function=NullFunc, take_log=True,
                convert_function=_get_convert("pres"),
                units=r"\rm{erg}\//\/\rm{cm}^{3}")
add_flash_field("pden", function=NullFunc, take_log=True,
                convert_function=_get_convert("pden"),
                units=r"\rm{g}/\rm{cm}^3")
add_flash_field("magx", function=NullFunc, take_log=False,
                convert_function=_get_convert("magx"),
                units = r"\mathrm{Gau\ss}")
add_flash_field("magy", function=NullFunc, take_log=False,
                convert_function=_get_convert("magy"),
                units = r"\mathrm{Gau\ss}")
add_flash_field("magz", function=NullFunc, take_log=False,
                convert_function=_get_convert("magz"),
                units = r"\mathrm{Gau\ss}")
add_flash_field("magp", function=NullFunc, take_log=True,
                convert_function=_get_convert("magp"),
                units = r"\rm{erg}\//\/\rm{cm}^{3}")
add_flash_field("divb", function=NullFunc, take_log=False,
                convert_function=_get_convert("divb"),
                units = r"\mathrm{Gau\ss}\/\rm{cm}")
add_flash_field("game", function=NullFunc, take_log=False,
                convert_function=_get_convert("game"),
                units=r"\rm{ratio\/of\/specific\/heats}")
add_flash_field("gamc", function=NullFunc, take_log=False,
                convert_function=_get_convert("gamc"),
                units=r"\rm{ratio\/of\/specific\/heats}")
add_flash_field("gpot", function=NullFunc, take_log=False,
                convert_function=_get_convert("gpot"),
                units=r"\rm{ergs\//\/g}")
add_flash_field("gpol", function=NullFunc, take_log=False,
                convert_function=_get_convert("gpol"),
                units = r"\rm{ergs\//\/g}")
add_flash_field("flam", function=NullFunc, take_log=False,
                convert_function=_get_convert("flam"))

for f,v in translation_dict.items():
    if v not in KnownFLASHFields:
        pfield = v.startswith("particle")
        add_flash_field(v, function=NullFunc, take_log=False,
                  validators = [ValidateDataField(v)],
                  particle_type = pfield)
    if f.endswith("_Fraction") :
        dname = "%s\/Fraction" % f.split("_")[0]
    else :
        dname = f                    
    ff = KnownFLASHFields[v]
    pfield = f.startswith("particle")
    add_field(f, TranslationFunc(v),
              take_log=KnownFLASHFields[v].take_log,
              units = ff._units, display_name=dname,
              particle_type = pfield)

def _convertParticleMassMsun(data):
    return 1.0/1.989e33
def _ParticleMassMsun(field, data):
    return data["ParticleMass"]
add_field("ParticleMassMsun",
          function=_ParticleMassMsun, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMassMsun,
          particle_convert_function=_ParticleMassMsun)

def _ThermalEnergy(fields, data) :
    try:
        return data["eint"]
    except:
        pass
    try:
        return data["Pressure"] / (data.pf["Gamma"] - 1.0) / data["Density"]
    except:
        pass
    if data.has_field_parameter("mu") :
        mu = data.get_field_parameter("mu")
    else:
        mu = 0.6
    return kboltz*data["Density"]*data["Temperature"]/(mu*mh) / (data.pf["Gamma"] - 1.0)
    
add_field("ThermalEnergy", function=_ThermalEnergy,
          units=r"\rm{ergs}/\rm{g}")

def _TotalEnergy(fields, data) :
    try:
        etot = data["ener"]
    except:
        etot = data["ThermalEnergy"] + 0.5 * (
            data["x-velocity"]**2.0 +
            data["y-velocity"]**2.0 +
            data["z-velocity"]**2.0)
    try:
        etot += data['magp']/data["Density"]
    except:
        pass
    return etot

add_field("TotalEnergy", function=_TotalEnergy,
          units=r"\rm{ergs}/\rm{g}")

def _GasEnergy(fields, data) :
    return data["ThermalEnergy"]

add_field("GasEnergy", function=_GasEnergy, 
          units=r"\rm{ergs}/\rm{g}")

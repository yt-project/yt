"""
Fields applicable only to Enzo

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

add_field = add_enzo_field

_speciesList = ["HI","HII","Electron",
               "HeI","HeII","HeIII",
               "H2I","H2II","HM",
               "DI","DII","HDI","Metal"]
def _SpeciesFraction(field, data):
    sp = field.name.split("_")[0] + "_Density"
    return data[sp]/data["Density"]
for species in _speciesList:
    add_field("%s_Fraction" % species,
             function=_SpeciesFraction,
             validators=ValidateDataField("%s_Density" % species))

def _Metallicity(field, data):
    return data["Metal_Fraction"] / 0.0204
add_field("Metallicity", units=r"Z_{\rm{Solar}}",
          function=_Metallicity,
          validators=ValidateDataField("Metal_Density"),
          projection_conversion="1")

def _ThermalEnergy(field, data):
    if data.pf["HydroMethod"] == 2:
        return data["Total_Energy"]
    elif data.pf["HydroMethod"] == 'orion':
        return data["Total_Energy"] - 0.5 * data["Density"] * (
                   data["x-velocity"]**2.0
                 + data["y-velocity"]**2.0
                 + data["z-velocity"]**2.0 )
    else:
        if data.pf["DualEnergyFormalism"]:
            return data["Gas_Energy"]
        else:
            return data["Total_Energy"] - 0.5*(
                   data["x-velocity"]**2.0
                 + data["y-velocity"]**2.0
                 + data["z-velocity"]**2.0 )
add_field("ThermalEnergy", function=_ThermalEnergy,
          units=r"\rm{ergs}/\rm{cm^3}")

def _NumberDensity(field, data):
    # We can assume that we at least have Density
    # We should actually be guaranteeing the presence of a .shape attribute,
    # but I am not currently implementing that
    fieldData = na.zeros(data["Density"].shape,
                         dtype = data["Density"].dtype)
    if data.pf["MultiSpecies"] == 0:
        if data.has_field_parameter("mu"):
            mu = data.get_field_parameter("mu")
        else:
            mu = 0.6
        fieldData += data["Density"] * mu
    if data.pf["MultiSpecies"] > 0:
        fieldData += data["HI_Density"] / 1.0
        fieldData += data["HII_Density"] / 1.0
        fieldData += data["HeI_Density"] / 4.0
        fieldData += data["HeII_Density"] / 4.0
        fieldData += data["HeIII_Density"] / 4.0
        fieldData += data["Electron_Density"] / 1.0
    if data.pf["MultiSpecies"] > 1:
        fieldData += data["HM_Density"] / 1.0
        fieldData += data["H2I_Density"] / 2.0
        fieldData += data["H2II_Density"] / 2.0
    if data.pf["MultiSpecies"] > 2:
        fieldData += data["DI_Density"] / 2.0
        fieldData += data["DII_Density"] / 2.0
        fieldData += data["HDI_Density"] / 3.0
    return fieldData
def _ConvertNumberDensity(data):
    return 1.0/mh
add_field("NumberDensity", units=r"\rm{cm}^{-3}",
          function=_NumberDensity,
          convert_function=_ConvertNumberDensity)

# Now we add all the fields that we want to control, but we give a null function
# This is every Enzo field we can think of.  This will be installation-dependent,
#if data.pf["HydroMethod"] == 'orion':
_default_fields = ["Density","Temperature","Gas_Energy","Total_Energy",
                   "x-velocity","y-velocity","z-velocity",
                   "x-momentum","y-momentum","z-momentum"]
# else:
#     _default_fields = ["Density","Temperature","Gas_Energy","Total_Energy",
#                        "x-velocity","y-velocity","z-velocity"]
_default_fields += [ "%s_Density" % sp for sp in _speciesList ]

for field in _default_fields:
    add_field(field, function=lambda a, b: None, take_log=True,
              validators=[ValidateDataField(field)], units=r"\rm{g}/\rm{cm}^3")
EnzoFieldInfo["x-velocity"].projection_conversion='1'
EnzoFieldInfo["y-velocity"].projection_conversion='1'
EnzoFieldInfo["z-velocity"].projection_conversion='1'

# Now we override

def _convertDensity(data):
    return data.convert("Density")
for field in ["Density"] + [ "%s_Density" % sp for sp in _speciesList ]:
    EnzoFieldInfo[field]._units = r"\rm{g}/\rm{cm}^3"
    EnzoFieldInfo[field]._projected_units = r"\rm{g}/\rm{cm}^2"
    EnzoFieldInfo[field]._convert_function=_convertDensity

add_field("Dark_Matter_Density", function=lambda a,b: None,
          convert_function=_convertDensity,
          validators=[ValidateDataField("Dark_Matter_Density"),
                      ValidateSpatial(0)],
          not_in_all = True)

def _convertEnergy(data):
    return data.convert("x-velocity")**2.0
EnzoFieldInfo["Gas_Energy"]._units = r"\rm{ergs}/\rm{g}"
EnzoFieldInfo["Gas_Energy"]._convert_function = _convertEnergy
EnzoFieldInfo["Total_Energy"]._units = r"\rm{ergs}/\rm{g}"
EnzoFieldInfo["Total_Energy"]._convert_function = _convertEnergy
EnzoFieldInfo["Temperature"]._units = r"\rm{K}"

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = EnzoFieldInfo["%s-velocity" % ax]
    f._units = r"\rm{cm}/\rm{s}"
    f._convert_function = _convertVelocity
    f.take_log = False


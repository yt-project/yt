"""
Athena-specific fields



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
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType, \
    NullFunc, \
    TranslationFunc
from yt.utilities.physical_constants import \
    kboltz,mh
import yt.data_objects.universal_fields

AthenaFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = AthenaFieldInfo.add_field

KnownAthenaFields = FieldInfoContainer()
add_athena_field = KnownAthenaFields.add_field

add_athena_field("density", function=NullFunc, take_log=False,
                 units=r"", projected_units =r"")

add_athena_field("pressure", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("total_energy", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("velocity_x", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("velocity_y", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("velocity_z", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("momentum_x", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("momentum_y", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("momentum_z", function=NullFunc, take_log=False,
                 units=r"")

add_athena_field("cell_centered_B_x", function=NullFunc, take_log=False,
                 units=r"", display_name=r"$\rm{cell\/centered\/B_x}$")

add_athena_field("cell_centered_B_y", function=NullFunc, take_log=False,
                 units=r"", display_name=r"$\rm{cell\/centered\/B_y}$")

add_athena_field("cell_centered_B_z", function=NullFunc, take_log=False,
                 units=r"", display_name=r"$\rm{cell\/centered\/B_z}$")

# In Athena, conservative or primitive variables may be written out.
# By default, yt concerns itself with primitive variables. The following
# field definitions allow for conversions to primitive variables in the
# case that the file contains the conservative ones.

def _convertDensity(data) :
    return data.convert("Density")
def _density(field, data) :
    return data["density"]
add_field("Density", function=_density, take_log=False,
          units=r"\rm{g}/\rm{cm}^3", projected_units=r"\rm{g}/\rm{cm}^2",
          convert_function=_convertDensity)

def _convertVelocity(data):
    return data.convert("x-velocity")
def _xvelocity(field, data):
    if "velocity_x" in data.pf.field_info:
        return data["velocity_x"]
    else:
        return data["momentum_x"]/data["density"]           
add_field("x-velocity", function=_xvelocity, take_log=False,
          units=r"\rm{cm}/\rm{s}", convert_function=_convertVelocity)
def _yvelocity(field, data):
    if "velocity_y" in data.pf.field_info:
        return data["velocity_y"]
    else:
        return data["momentum_y"]/data["density"]
add_field("y-velocity", function=_yvelocity, take_log=False,
          units=r"\rm{cm}/\rm{s}", convert_function=_convertVelocity)
def _zvelocity(field, data):
    if "velocity_z" in data.pf.field_info:
        return data["velocity_z"]
    else:
        return data["momentum_z"]/data["density"]
add_field("z-velocity", function=_zvelocity, take_log=False,
          units=r"\rm{cm}/\rm{s}", convert_function=_convertVelocity)

def _convertEnergy(data) :
    return data.convert("x-velocity")**2
def _gasenergy(field, data) :
    if "pressure" in data.pf.field_info:
        return data["pressure"]/(data.pf["Gamma"]-1.0)/data["density"]
    else:
        eint = data["total_energy"] - 0.5*(data["momentum_x"]**2 +
                                           data["momentum_y"]**2 +
                                           data["momentum_z"]**2)/data["density"]
        if "cell_centered_B_x" in data.pf.field_info:
            eint -= 0.5*(data["cell_centered_B_x"]**2 +
                         data["cell_centered_B_y"]**2 +
                         data["cell_centered_B_z"]**2)
        return eint/data["density"]
add_field("Gas_Energy", function=_gasenergy, take_log=False,
          convert_function=_convertEnergy, units=r"\rm{erg}/\rm{g}")

def _convertPressure(data) :
    return data.convert("Density")*data.convert("x-velocity")**2
def _pressure(field, data) :
    if "pressure" in data.pf.field_info:
        return data["pressure"]
    else:
        eint = data["total_energy"] - 0.5*(data["momentum_x"]**2 +
                                           data["momentum_y"]**2 +
                                           data["momentum_z"]**2)/data["density"]
        if "cell_centered_B_x" in data.pf.field_info:
            eint -= 0.5*(data["cell_centered_B_x"]**2 +
                         data["cell_centered_B_y"]**2 +
                         data["cell_centered_B_z"]**2)
        return eint*(data.pf["Gamma"]-1.0)
add_field("Pressure", function=_pressure, take_log=False,
          convert_function=_convertPressure, units=r"\rm{erg}/\rm{cm}^3",
          projected_units=r"\rm{erg}/\rm{cm}^2")

def _temperature(field, data):
    if data.has_field_parameter("mu"):
        mu = data.get_field_parameter("mu")
    else:
        mu = 0.6
    return mu*mh*data["Pressure"]/data["Density"]/kboltz
add_field("Temperature", function=_temperature, take_log=False,
          units=r"\rm{K}")

def _convertBfield(data):
        return np.sqrt(4*np.pi*data.convert("Density")*data.convert("x-velocity")**2)
def _Bx(field, data):
    return data['cell_centered_B_x']
add_field("Bx", function=_Bx, take_log=False,
          units=r"\rm{Gauss}", display_name=r"B_x",
          convert_function=_convertBfield)
def _By(field, data):
    return data['cell_centered_B_y']
add_field("By", function=_By, take_log=False,
          units=r"\rm{Gauss}", display_name=r"B_y",
          convert_function=_convertBfield)
def _Bz(field, data):
    return data['cell_centered_B_z']
add_field("Bz", function=_Bz, take_log=False,
          units=r"\rm{Gauss}", display_name=r"B_z",
          convert_function=_convertBfield)



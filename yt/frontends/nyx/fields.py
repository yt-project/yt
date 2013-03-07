"""
Field specifications for Nyx

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Casey W. Stark, J. S. Oishi, Matthew Turk.  All Rights
  Reserved.

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

import yt.data_objects.universal_fields

from yt.data_objects.field_info_container import FieldInfoContainer, \
    NullFunc, TranslationFunc, FieldInfo, \
    ValidateParameter, ValidateDataField, ValidateProperty, ValidateSpatial, \
    ValidateGridType
from yt.utilities.physical_constants import mh, kboltz

NyxFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = NyxFieldInfo.add_field

KnownNyxFields = FieldInfoContainer()
add_nyx_field = KnownNyxFields.add_field

# Density
add_nyx_field("density", function=lambda a, b: None, take_log=True,
          validators=[ValidateDataField("density")],
          units="g/cm**3",
          projected_units =r"\rm{g} / \rm{cm}^2")
KnownNyxFields["density"]._projected_units =r"\rm{g} / \rm{cm}^2"

add_field("Density", function=TranslationFunc("density"), take_log=True,
          units=r"\rm{g} / \rm{cm}^3",
          projected_units =r"\rm{g} / \rm{cm}^2")

# Particle mass in units of $ M_{\odot}
def _convertParticleMassMsun(data):
    return (1/1.989e33)
def _particle_mass_m_sun(field, data):
    return data["particle_mass"]
add_field("ParticleMassMsun", function=_particle_mass_m_sun,
          validators=[ValidateSpatial(0), ValidateDataField("particle_mass")],
          particle_type=True, convert_function=_convertParticleMassMsun,
          take_log=True, units=r"\rm{M_{\odot}}")

add_nyx_field("Dark_Matter_Density", function=TranslationFunc("particle_mass_density"),
          take_log=True,
          units=r"\rm{g} / \rm{cm}^3",particle_type=True,
          projected_units =r"\rm{g} / \rm{cm}^2")


# Energy Density
# @todo: ``energy_density``
add_nyx_field("total_energy", function=lambda a, b: None, take_log=True,
          validators=[ValidateDataField("total_energy")],
          units=r"\rm{M_{\odot}} (\rm{km} / \rm{s})^2")

# Momentum in each dimension.
# @todo: ``momentum_x``
add_nyx_field("x-momentum", function=lambda a, b: None, take_log=False,
          validators=[ValidateDataField("x-momentum")],
          units=r"\rm{M_{\odot}} \rm{km} / \rm{s}")
add_nyx_field("y-momentum", function=lambda a, b: None, take_log=False,
          validators=[ValidateDataField("y-momentum")],
          units=r"\rm{M_{\odot}} \rm{km} / \rm{s}")
add_nyx_field("z-momentum", function=lambda a, b: None, take_log=False,
          validators=[ValidateDataField("z-momentum")],
          units=r"\rm{M_{\odot}} \rm{km} / \rm{s}")

### Now derived fields

# Velocity fields in each dimension
# @todo: ``velocity_x``
def _x_velocity(field, data):
    """ Generate x-velocity from x-momentum and density. """
    return data["x-momentum"] / data["density"]
add_field("x-velocity", function=_x_velocity, take_log=False,
          units=r"\rm{km} / \rm{s}")

def _y_velocity(field, data):
    """ Generate y-velocity from y-momentum and density. """
    return data["y-momentum"] / data["density"]
add_field("y-velocity", function=_y_velocity, take_log=False,
          units=r"\rm{km} / \rm{s}")

def _z_velocity(field, data):
    """ Generate z-velocity from z-momentum and density. """
    return data["z-momentum"] / data["density"]
add_field("z-velocity", function=_z_velocity, take_log=False,
          units=r"\rm{km} / \rm{s}")

# The gas **thermal** energy.
# @todo: should be called ``gas_energy`` whether it is data or derived
def _thermal_energy(field, data):
    """
    Generate thermal (gas energy). Dual Energy Formalism was implemented by
    Stella, but this isn't how it's called, so I'll leave that commented out for
    now.

    """
    #if data.pf["DualEnergyFormalism"]:
    #    return data["Gas_Energy"]
    #else:
    return data["Total_Energy"] - 0.5 * data["density"] * (
                                          data["x-velocity"]**2.0
                                        + data["y-velocity"]**2.0
                                        + data["z-velocity"]**2.0 )
add_field("ThermalEnergy", function=_thermal_energy,
          units=r"\rm{M_{\odot}} (\rm{km} / \rm{s})^2")

# Gas pressure
# @todo: eventually figure out a way to detect when using radiation and change
#        this.
def _pressure(field, data):
    """
    Computed using

    $$ pressure = (\gamma - 1.0) * e$$

    where e is thermal energy density. Note that this will need to be modified
    when radiation is accounted for.

    """
    return (data.pf.gamma - 1.0) * data["ThermalEnergy"]
add_field("Pressure", function=_pressure,
          units=r"\rm{M_{\odot}} (\rm{km} / \rm{s})^2 / \rm{Mpc}^3")

# Gas temperature
def _temperature(field, data):
    return ((data.pf.gamma - 1.0) * data.pf["mu"] * mh *
            data["ThermalEnergy"] / (kboltz * data["Density"]))
add_field("Temperature", function=_temperature, take_log=False,
          units=r"\rm{Kelvin}")


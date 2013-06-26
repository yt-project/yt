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

#
# Constants
#

# BEWARE hardcoded, uniform gamma value
nyx_gamma = 5.0 / 3.0

# Density
add_nyx_field("density", function=lambda a, b: None, take_log=True,
              validators=[ValidateDataField("density")], units="g/cm**3")

add_field("Density", function=TranslationFunc("density"), take_log=True,
          units="g/cm**3")

# Particle mass in units of $ M_{\odot}
def _convertParticleMassMsun(data):
    return (1/1.989e33)
def _particle_mass_m_sun(field, data):
    return data["particle_mass"]
add_field("ParticleMassMsun", function=_particle_mass_m_sun,
          validators=[ValidateSpatial(0), ValidateDataField("particle_mass")],
          particle_type=True, convert_function=_convertParticleMassMsun,
          take_log=True, units="Msun")

add_nyx_field("Dark_Matter_Density",
              function=TranslationFunc("particle_mass_density"),
              take_log=True, particle_type=True, units="g/cm**3")


add_nyx_field("energy_density", function=lambda a, b: None, take_log=True,
              validators=[ValidateDataField("total_energy")],
              units="Msun * (km / s)**2")

add_nyx_field("momentum_x", function=lambda a, b: None, take_log=False,
              validators=[ValidateDataField("x-momentum")], units="Msun*km/s")
add_nyx_field("momentum_y", function=lambda a, b: None, take_log=False,
              validators=[ValidateDataField("y-momentum")], units="Msun*km/s")
add_nyx_field("momentum_z", function=lambda a, b: None, take_log=False,
              validators=[ValidateDataField("z-momentum")], units="Msun*km/s")

### Now derived fields

# Velocity fields in each dimension
# @todo: ``velocity_x``
def _velocity_x(field, data):
    """ Generate x-velocity from x-momentum and density. """
    return data["momentum_x"] / data["density"]
add_field("velocity_x", function=_velocity_x, take_log=False, units="km/s")

def _velocity_y(field, data):
    """ Generate y-velocity from y-momentum and density. """
    return data["momentum_y"] / data["density"]
add_field("velocity_y", function=_velocity_y, take_log=False, units="km/s")

def _velocity_z(field, data):
    """ Generate z-velocity from z-momentum and density. """
    return data["momentum_z"] / data["density"]
add_field("velocity_z", function=_velocity_z, take_log=False, units="km/s")

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
    return ( data["total_energy"]
             - 0.5 * data["density"] * (   data["velocity_x"]**2.0
                                         + data["velocity_y"]**2.0
                                         + data["velocity_z"]**2.0 ) )
add_field("thermal_energy", function=_thermal_energy, units="Msun*(km/s)**2")

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
    return (nyx_gamma - 1.0) * data["ThermalEnergy"]

add_field("pressure", function=_pressure, units="Msun*(km/s)**2/Mpc**3")

# Gas temperature
def _temperature(field, data):
    return ( (gamma - 1.0) * data.pf["mu"] * mh *
             data["thermal_energy"] / (kboltz * data["density"]) )
add_field("temperature", function=_temperature, take_log=False, units="K")

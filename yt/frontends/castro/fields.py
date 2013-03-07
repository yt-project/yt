"""
Castro-specific fields

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: UC Berkeley
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2010 J. S. Oishi, Matthew Turk.  All Rights Reserved.

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
    TranslationFunc, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
from yt.utilities.physical_constants import mh, kboltz

translation_dict = {
    "x-velocity": "xvel",
    "y-velocity": "yvel",
    "z-velocity": "zvel",
    "Density": "density",
    "Total_Energy": "eden",
    "Temperature": "temperature",
    "x-momentum": "xmom",
    "y-momentum": "ymom",
    "z-momentum": "zmom"
}

# Setup containers for fields possibly in the output files
KnownCastroFields = FieldInfoContainer()
add_castro_field = KnownCastroFields.add_field

# and always derived ones
CastroFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = CastroFieldInfo.add_field

# Start adding fields
add_castro_field("density", function=NullFunc, take_log=True,
                 units="g/cm**3")

add_castro_field("eden", function=NullFunc, take_log=True,
                 validators = [ValidateDataField("eden")],
                 units="erg/cm**3")

add_castro_field("xmom", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("xmom")],
                 units="g/cm**2/s")

add_castro_field("ymom", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("ymom")],
                 units="g/cm**2/s")

add_castro_field("zmom", function=NullFunc, take_log=False,
                 validators = [ValidateDataField("zmom")],
                 units="g/cm**2/s")

# Now populate derived fields
for mine, theirs in translation_dict.items():
    if KnownCastroFields.has_key(theirs):
        add_field(theirs, function=TranslationFunc(mine),
                  take_log=KnownCastroFields[theirs].take_log)

# Now fallbacks, in case these fields are not output
def _xVelocity(field, data):
    """ Generate x-velocity from x-momentum and density. """
    return data["xmom"] / data["density"]

add_field("x-velocity", function=_xVelocity, take_log=False,
          units='cm/s')

def _yVelocity(field, data):
    """ Generate y-velocity from y-momentum and density. """
    return data["ymom"] / data["density"]

add_field("y-velocity", function=_yVelocity, take_log=False,
          units='cm/s')

def _zVelocity(field, data):
    """ Generate z-velocity from z-momentum and density. """
    return data["zmom"] / data["density"]

add_field("z-velocity", function=_zVelocity, take_log=False,
          units='cm/s')

def _ThermalEnergy(field, data):
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

add_field("ThermalEnergy", function=_ThermalEnergy,
          units="erg/cm**3")

def _Pressure(field, data):
    """
    M{(Gamma-1.0)*e, where e is thermal energy density
    
    NB: this will need to be modified for radiation

    """
    return (data.pf.gamma - 1.0) * data["ThermalEnergy"]

add_field("Pressure", function=_Pressure, units="dyne/cm**2")

def _Temperature(field, data):
    return ((data.pf.gamma - 1.0) * data.pf["mu"] * mh *
            data["ThermalEnergy"] / (kboltz * data["Density"]))

add_field("Temperature", function=_Temperature, units="K",
          take_log=False)

def _convertParticleMassMsun(data):
    return 1.0 / 1.989e33
def _ParticleMassMsun(field, data):
    return data["particle_mass"]

add_field("ParticleMassMsun",
          function=_ParticleMassMsun, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMassMsun,
          particle_convert_function=_ParticleMassMsun)

# Fundamental fields that are usually/always output:
#   density
#   xmom
#   ymom
#   zmom
#   rho_E
#   rho_e
#   Temp
#
# "Derived" fields that are sometimes output:
#   x_velocity
#   y_velocity
#   z_velocity
#   magvel
#   grav_x
#   grav_y
#   grav_z
#   maggrav
#   magvort
#   pressure
#   entropy
#   divu
#   eint_e (e as derived from the "rho e" variable)
#   eint_E (e as derived from the "rho E" variable)

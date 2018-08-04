"""
Cholla-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.utilities.physical_constants import \
    kboltz, mh

erg_units = "code_mass * (code_length*code_time**2)"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
vel_units = "code_length / code_time"

def velocity_field(comp):
    def _velocity(field, data):
        return data["cholla", "momentum_%s" % comp]/data["cholla","density"]
    return _velocity

class ChollaFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", (rho_units, ["density"], None)),
        ("Energy", (erg_units, ["total_energy"], None)),
        ("GasEnergy", (erg_units, ["gas_energy"], None)),
        ("momentum_x", (mom_units, ["momentum_x"], None)),
        ("momentum_y", (mom_units, ["momentum_y"], None)),
        ("momentum_z", (mom_units, ["momentum_z"], None)),
    )

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system
        # Add velocity fields
        for comp in "xyz":
            mom_field = ("cholla", "momentum_%s" % comp)
            if mom_field in self.field_list:
                self.add_output_field(mom_field, sampling_type="cell",
                                      units=mom_units)
                self.add_field(("gas","velocity_%s" % comp), sampling_type="cell",
                               function=velocity_field(comp), units = unit_system["velocity"])
        # Add pressure, energy, and temperature fields
        self.add_output_field(('cholla', 'gas_energy'), sampling_type='cell',
                              units=erg_units)
        self.alias(('gas', 'thermal_energy'), ('cholla', 'gas_energy'),
                   units=unit_system['energy'])
        self.add_output_field(('cholla', 'energy'), sampling_type='cell',
                              units=erg_units)
        self.alias(('gas', 'total_energy'), ('cholla', 'energy'),
                   units=unit_system['energy'])
        # Add pressure field            
        def _pressure(field, data):
            return data[('cholla', 'gas_energy')] * (data.ds.gamma-1)
        self.add_field(("gas","pressure"), sampling_type="cell",  function=_pressure,
                       units=unit_system["pressure"])
        # Add temperature field
        def _temperature(field, data):
            if data.has_field_parameter("mu"):
                mu = data.get_field_parameter("mu")
            else:
                mu = 0.6
            return mu*mh*data["gas","pressure"]/data["gas","density"]/kboltz
        self.add_field(("gas","temperature"), sampling_type="cell",  function=_temperature,
                       units=unit_system["temperature"])

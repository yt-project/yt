"""
GAMER-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import mh, boltzmann_constant_cgs

b_units   = "code_magnetic"
pre_units = "code_mass / (code_length*code_time**2)"
erg_units = "code_mass / (code_length*code_time**2)"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_length**2*code_time)"
vel_units = "code_velocity"
pot_units = "code_length**2/code_time**2"

psi_units = "code_mass**(1/2) / code_length**(3/2)"


class GAMERFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # hydro fields on disk (GAMER outputs conservative variables)
        ( "Dens", (rho_units, ["density"],                 r"\rho") ),
        ( "MomX", (mom_units, ["momentum_x"],              None   ) ),
        ( "MomY", (mom_units, ["momentum_y"],              None   ) ),
        ( "MomZ", (mom_units, ["momentum_z"],              None   ) ),
        ( "Engy", (erg_units, ["total_energy_per_volume"], None   ) ),
        ( "Pote", (pot_units, ["gravitational_potential"], None   ) ),

        # psiDM fields on disk
        ( "Real", (psi_units, ["psidm_real_part"],         None   ) ),
        ( "Imag", (psi_units, ["psidm_imaginary_part"],    None   ) ),
    )

    known_particle_fields = (
    )

    def __init__(self, ds, field_list):
        super(GAMERFieldInfo, self).__init__(ds, field_list)

    # add primitive and other derived variables
    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system

        # velocity
        def velocity_xyz(v):
            def _velocity(field, data):
                return data["gas", "momentum_%s"%v] / data["gas","density"]
            return _velocity
        for v in "xyz":
            self.add_field( ("gas","velocity_%s"%v), function = velocity_xyz(v),
                            units = unit_system["velocity"] )

        # ============================================================================
        # note that yt internal fields assume
        #    [thermal_energy]          = [energy per mass]
        #    [kinetic_energy]          = [energy per volume]
        # and we further adopt
        #    [total_energy]            = [energy per mass]
        #    [total_energy_per_volume] = [energy per volume]
        # ============================================================================

        # kinetic energy per volume
        def ek(data):
            return 0.5*( data["gamer","MomX"]**2 +
                         data["gamer","MomY"]**2 +
                         data["gamer","MomZ"]**2 ) / data["gamer","Dens"]

        # thermal energy per volume
        def et(data):
            return data["gamer","Engy"] - ek(data)

        # thermal energy per mass (i.e., specific)
        def _thermal_energy(field, data):
            return et(data) / data["gamer","Dens"]
        self.add_field( ("gas","thermal_energy"), function = _thermal_energy,
                        units = unit_system["specific_energy"] )

        # total energy per mass
        def _total_energy(field, data):
            return data["gamer","Engy"] / data["gamer","Dens"]
        self.add_field( ("gas","total_energy"), function = _total_energy,
                        units = unit_system["specific_energy"] )

        # pressure
        def _pressure(field, data):
            return et(data)*(data.ds.gamma-1.0)
        self.add_field( ("gas","pressure"), function = _pressure,
                        units = unit_system["pressure"] )

        # temperature
        def _temperature(field, data):
            return data.ds.mu*mh*data["gas","pressure"] / \
                   (data["gas","density"]*boltzmann_constant_cgs)
        self.add_field( ("gas","temperature"), function = _temperature,
                        units = unit_system["temperature"] )

    def setup_particle_fields(self, ptype):
        # This will get called for every particle type.
        pass

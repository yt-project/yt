"""
AdaptaHOP-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

m_units = "1e11 * Msun"
r_units = "Mpc"
v_units = "km / s"
l_units = "1e11 * Msun * Mpc * km / s"
e_units = "1e11 * Msun * km**2 / s**2"
dens_units = "1e11 * Msun / Mpc**3"

class AdaptaHOPFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("particle_identifier", ("", [], "Halo Identity")),
        ("raw_position_x", (r_units, [], None)),
        ("raw_position_y", (r_units, [], None)),
        ("raw_position_z", (r_units, [], None)),
        ("r", (r_units, [], "Halo radius")),
        ("a", (r_units, [], "Halo semi-major axis")),
        ("b", (r_units, [], "Halo semi-medium axis")),
        ("c", (r_units, [], "Halo semi-minor axis")),
        ("particle_velocity_x", (v_units, [], "Halo velocity x")),
        ("particle_velocity_y", (v_units, [], "Halo velocity y")),
        ("particle_velocity_z", (v_units, [], "Halo velocity z")),
        ("particle_angular_momentum_x", (l_units, [], "Halo Angular Momentum x")),
        ("particle_angular_momentum_y", (l_units, [], "Halo Angular Momentum y")),
        ("particle_angular_momentum_z", (l_units, [], "Halo Angular Momentum z")),
        ("particle_mass", (m_units, [], "Halo mass")),
        ("ek", (e_units, [], "Halo Kinetic Energy")),
        ("ep", (e_units, [], "Halo Gravitational Energy")),
        ("ek", (e_units, [], "Halo Total Energy")),
        ("spin", ("", [], "Halo Spin")),
        # Virial parameters
        ("virial_radius", (r_units, [], "Halo Virial Radius")),
        ("virial_mass", (m_units, [], "Halo Virial Mass")),
        ("virial_temperature", ("K", [], "Halo Virial Temperature")),
        ("virial_velocity", (v_units, [], "Halo Virial Velocity")),
        # NFW parameters
        ("rho0", (dens_units, [], "Halo NFW Density")),
        ("R_c", (dens_units, [], "Halo NFW Scale Radius")),
    )

    def setup_particle_fields(self, ptype):
        super(AdaptaHOPFieldInfo, self).setup_particle_fields(ptype)

        # Add particle position
        def generate_pos_field(d):
            shift = self.ds.domain_width[0] / 2
            def closure(field, data):
                return data["halos", "raw_position_%s" % d] + shift
            return closure

        for k in 'xyz':
            fun = generate_pos_field(k)
            self.add_field(("halos", "particle_position_%s" % k), sampling_type="particle", function=fun,
                           units='Mpc')

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

m_units = "1e11 * Msun"       # Msun / h
r_units = "Mpc"               # Mpc / h comoving
v_units = "km / s"            # km /s phys, peculiar

class AdaptaHOPFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("particle_identifier", ("", [], None)),
        ("raw_position_x", (r_units, [], None)),
        ("raw_position_y", (r_units, [], None)),
        ("raw_position_z", (r_units, [], None)),
        ("particle_velocity_x", (v_units, [], None)),
        ("particle_velocity_y", (v_units, [], None)),
        ("particle_velocity_z", (v_units, [], None)),
        ("particle_corevel_x", (v_units, [], None)),
        ("particle_corevel_y", (v_units, [], None)),
        ("particle_corevel_z", (v_units, [], None)),
        ("particle_bulkvel_x", (v_units, [], None)),
        ("particle_bulkvel_y", (v_units, [], None)),
        ("particle_bulkvel_z", (v_units, [], None)),
        ("particle_mass", (m_units, [], "Mass")),
        ("virial_radius", (r_units, [], "Radius")),
        ("child_r", (r_units, [], None)),
        ("vmax_r", (v_units, [], None)),
        # These fields I don't have good definitions for yet.
        ('mgrav', ("", [], None)),
        ('vmax', (v_units, [], "V_{max}")),
        ('rvmax', (v_units, [], None)),
        ('rs', (r_units, [], "R_s")),
        ('klypin_rs', (r_units, [], "Klypin R_s")),
        ('vrms', (v_units, [], "V_{rms}")),
        ('Jx', ("", [], "J_x")),
        ('Jy', ("", [], "J_y")),
        ('Jz', ("", [], "J_z")),
        ('energy', ("", [], None)),
        ('spin', ("", [], "Spin Parameter")),
        ('alt_m1', (m_units, [], None)),
        ('alt_m2', (m_units, [], None)),
        ('alt_m3', (m_units, [], None)),
        ('alt_m4', (m_units, [], None)),
        ('Xoff', ("", [], None)),
        ('Voff', ("", [], None)),
        ('b_to_a', ("", [], "Ellipsoidal b to a")),
        ('c_to_a', ("", [], "Ellipsoidal c to a")),
        ('Ax', ("", [], "A_x")),
        ('Ay', ("", [], "A_y")),
        ('Az', ("", [], "A_z")),
        ('b_to_a2', ("", [], None)),
        ('c_to_a2', ("", [], None)),
        ('A2x', ("", [], "A2_x")),
        ('A2y', ("", [], "A2_y")),
        ('A2z', ("", [], "A2_z")),
        ('bullock_spin', ("", [], "Bullock Spin Parameter")),
        ('kin_to_pot', ("", [], "Kinetic to Potential")),
        ('m_pe_b', ("", [], None)),
        ('m_pe_d', ("", [], None)),
        ('num_p', ("", [], "Number of Particles")),
        ('num_child_particles', ("", [], "Number of Child Particles")),
        ('p_start', ("", [], None)),
        ('desc', ("", [], None)),
        ('flags', ("", [], None)),
        ('n_core', ("", [], None)),
        ('min_pos_err', ("", [], None)),
        ('min_vel_err', ("", [], None)),
        ('min_bulkvel_err', ("", [], None)),
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

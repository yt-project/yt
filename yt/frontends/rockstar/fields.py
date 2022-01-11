from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer

m_units = "Msun / h"  # Msun / h
p_units = "Mpccm / h"  # Mpc / h comoving
v_units = "km / s"  # km /s phys, peculiar
r_units = "kpccm / h"  # kpc / h comoving


class RockstarFieldInfo(FieldInfoContainer):

    known_particle_fields: KnownFieldsT = (
        ("particle_identifier", ("", [], None)),
        ("particle_position_x", (p_units, [], None)),
        ("particle_position_y", (p_units, [], None)),
        ("particle_position_z", (p_units, [], None)),
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
        ("mgrav", ("", [], None)),
        ("vmax", (v_units, [], "V_{max}")),
        ("rvmax", (v_units, [], None)),
        ("rs", (r_units, [], "R_s")),
        ("klypin_rs", (r_units, [], "Klypin R_s")),
        ("vrms", (v_units, [], "V_{rms}")),
        ("Jx", ("", [], "J_x")),
        ("Jy", ("", [], "J_y")),
        ("Jz", ("", [], "J_z")),
        ("energy", ("", [], None)),
        ("spin", ("", [], "Spin Parameter")),
        ("alt_m1", (m_units, [], None)),
        ("alt_m2", (m_units, [], None)),
        ("alt_m3", (m_units, [], None)),
        ("alt_m4", (m_units, [], None)),
        ("Xoff", ("", [], None)),
        ("Voff", ("", [], None)),
        ("b_to_a", ("", [], "Ellipsoidal b to a")),
        ("c_to_a", ("", [], "Ellipsoidal c to a")),
        ("Ax", ("", [], "A_x")),
        ("Ay", ("", [], "A_y")),
        ("Az", ("", [], "A_z")),
        ("b_to_a2", ("", [], None)),
        ("c_to_a2", ("", [], None)),
        ("A2x", ("", [], "A2_x")),
        ("A2y", ("", [], "A2_y")),
        ("A2z", ("", [], "A2_z")),
        ("bullock_spin", ("", [], "Bullock Spin Parameter")),
        ("kin_to_pot", ("", [], "Kinetic to Potential")),
        ("m_pe_b", ("", [], None)),
        ("m_pe_d", ("", [], None)),
        ("num_p", ("", [], "Number of Particles")),
        ("num_child_particles", ("", [], "Number of Child Particles")),
        ("p_start", ("", [], None)),
        ("desc", ("", [], None)),
        ("flags", ("", [], None)),
        ("n_core", ("", [], None)),
        ("min_pos_err", ("", [], None)),
        ("min_vel_err", ("", [], None)),
        ("min_bulkvel_err", ("", [], None)),
    )

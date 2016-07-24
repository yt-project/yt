"""
OWLSSubfind-specific fields




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

m_units = "code_mass"
mdot_units = "code_mass / code_time"
p_units = "Mpccm/h"
v_units = "1e5 * cmcm / s"

class OWLSSubfindFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("CenterOfMass_0", (p_units, ["particle_position_x"], None)),
        ("CenterOfMass_1", (p_units, ["particle_position_y"], None)),
        ("CenterOfMass_2", (p_units, ["particle_position_z"], None)),
        ("CenterOfMassVelocity_0", (v_units, ["particle_velocity_x"], None)),
        ("CenterOfMassVelocity_1", (v_units, ["particle_velocity_y"], None)),
        ("CenterOfMassVelocity_2", (v_units, ["particle_velocity_z"], None)),
        ("Mass", (m_units, ["particle_mass"], None)),
        ("Halo_M_Crit200", (m_units, ["Virial Mass"], None)),
        ("Halo_M_Crit2500", (m_units, [], None)),
        ("Halo_M_Crit500", (m_units, [], None)),
        ("Halo_M_Mean200", (m_units, [], None)),
        ("Halo_M_Mean2500", (m_units, [], None)),
        ("Halo_M_Mean500", (m_units, [], None)),
        ("Halo_M_TopHat200", (m_units, [], None)),
        ("Halo_R_Crit200", (p_units, ["Virial Radius"], None)),
        ("Halo_R_Crit2500", (p_units, [], None)),
        ("Halo_R_Crit500", (p_units, [], None)),
        ("Halo_R_Mean200", (p_units, [], None)),
        ("Halo_R_Mean2500", (p_units, [], None)),
        ("Halo_R_Mean500", (p_units, [], None)),
        ("Halo_R_TopHat200", (p_units, [], None)),
        ("BH_Mass", (m_units, [], None)),
        ("Stars/Mass", (m_units, [], None)),
        ("BH_Mdot", (mdot_units, [], None)),
        ("StarFormationRate", (mdot_units, [], None)),
)

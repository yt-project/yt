"""
GadgetFOF-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

m_units = "code_mass"
p_units = "code_length"
v_units = "code_velocity"

class GadgetFOFFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("GroupPos_0", (p_units, ["Group", "particle_position_x"], None)),
        ("GroupPos_1", (p_units, ["Group", "particle_position_y"], None)),
        ("GroupPos_2", (p_units, ["Group", "particle_position_z"], None)),
        ("GroupVel_0", (v_units, ["Group", "particle_velocity_x"], None)),
        ("GroupVel_1", (v_units, ["Group", "particle_velocity_y"], None)),
        ("GroupVel_2", (v_units, ["Group", "particle_velocity_z"], None)),
        ("GroupMass",  (m_units, ["Group", "particle_mass"], None)),
        ("GroupLen",   ("",      ["Group", "particle_number"], None)),
        ("SubhaloPos_0", (p_units, ["Subhalo", "particle_position_x"], None)),
        ("SubhaloPos_1", (p_units, ["Subhalo", "particle_position_y"], None)),
        ("SubhaloPos_2", (p_units, ["Subhalo", "particle_position_z"], None)),
        ("SubhaloVel_0", (v_units, ["Subhalo", "particle_velocity_x"], None)),
        ("SubhaloVel_1", (v_units, ["Subhalo", "particle_velocity_y"], None)),
        ("SubhaloVel_2", (v_units, ["Subhalo", "particle_velocity_z"], None)),
        ("SubhaloMass",  (m_units, ["Subhalo", "particle_mass"], None)),
        ("SubhaloLen",   ("",      ["Subhalo", "particle_number"], None)),
    )

    # these are extra fields to be created for the "all" particle type
    extra_union_fields = (
        (p_units, "particle_position_x"),
        (p_units, "particle_position_y"),
        (p_units, "particle_position_z"),
        (v_units, "particle_velocity_x"),
        (v_units, "particle_velocity_y"),
        (v_units, "particle_velocity_z"),
        (m_units, "particle_mass"),
        ("",      "particle_number"),
        ("",      "particle_ones"),
    )

class GadgetFOFHaloFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("ID", ("", ["Group", "particle_identifier"], None)),
    )

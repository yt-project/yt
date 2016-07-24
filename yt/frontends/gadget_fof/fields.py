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

_pnums = 6
_type_fields = \
  tuple(("%s%sType_%d" % (ptype, field, pnum), (units, [], None))
        for pnum in range(_pnums)
        for field, units in (("Mass", m_units), ("Len", p_units))
        for ptype in ("Group", "Subhalo"))
_sub_type_fields = \
  tuple(("Subhalo%sType_%d" % (field, pnum), (units, [], None))
        for pnum in range(_pnums)
        for field, units in (("HalfmassRad", p_units),
                             ("MassInHalfRad", m_units),
                             ("MassInMaxRad", m_units),
                             ("MassInRad", m_units)))

_particle_fields = (
    ("GroupPos_0",           (p_units, ["Group", "particle_position_x"], None)),
    ("GroupPos_1",           (p_units, ["Group", "particle_position_y"], None)),
    ("GroupPos_2",           (p_units, ["Group", "particle_position_z"], None)),
    ("GroupVel_0",           (v_units, ["Group", "particle_velocity_x"], None)),
    ("GroupVel_1",           (v_units, ["Group", "particle_velocity_y"], None)),
    ("GroupVel_2",           (v_units, ["Group", "particle_velocity_z"], None)),
    ("GroupMass",            (m_units, ["Group", "particle_mass"], None)),
    ("GroupLen",             ("",      ["Group", "particle_number"], None)),
    ("GroupNsubs",           ("",      ["Group", "subhalo_number"], None)),
    ("GroupFirstSub",        ("",      [], None)),
    ("Group_M_Crit200",      (m_units, [], None)),
    ("Group_M_Crit500",      (m_units, [], None)),
    ("Group_M_Mean200",      (m_units, [], None)),
    ("Group_M_TopHat200",    (m_units, [], None)),
    ("Group_R_Crit200",      (p_units, [], None)),
    ("Group_R_Crit500",      (p_units, [], None)),
    ("Group_R_Mean200",      (p_units, [], None)),
    ("Group_R_TopHat200",    (p_units, [], None)),
    ("SubhaloPos_0",         (p_units, ["Subhalo", "particle_position_x"], None)),
    ("SubhaloPos_1",         (p_units, ["Subhalo", "particle_position_y"], None)),
    ("SubhaloPos_2",         (p_units, ["Subhalo", "particle_position_z"], None)),
    ("SubhaloVel_0",         (v_units, ["Subhalo", "particle_velocity_x"], None)),
    ("SubhaloVel_1",         (v_units, ["Subhalo", "particle_velocity_y"], None)),
    ("SubhaloVel_2",         (v_units, ["Subhalo", "particle_velocity_z"], None)),
    ("SubhaloMass",          (m_units, ["Subhalo", "particle_mass"], None)),
    ("SubhaloLen",           ("",      ["Subhalo", "particle_number"], None)),
    ("SubhaloCM_0",          (p_units, ["Subhalo", "center_of_mass_x"], None)),
    ("SubhaloCM_1",          (p_units, ["Subhalo", "center_of_mass_y"], None)),
    ("SubhaloCM_2",          (p_units, ["Subhalo", "center_of_mass_z"], None)),
    ("SubhaloSpin_0",        ("",      ["Subhalo", "spin_x"], None)),
    ("SubhaloSpin_1",        ("",      ["Subhalo", "spin_y"], None)),
    ("SubhaloSpin_2",        ("",      ["Subhalo", "spin_z"], None)),
    ("SubhaloGrNr",          ("",      ["Subhalo", "group_identifier"], None)),
    ("SubhaloHalfmassRad",   (p_units, [], None)),
    ("SubhaloIDMostbound",   ("",      [], None)),
    ("SubhaloMassInHalfRad", (m_units, [], None)),
    ("SubhaloMassInMaxRad",  (m_units, [], None)),
    ("SubhaloMassInRad",     (m_units, [], None)),
    ("SubhaloParent",        ("",      [], None)),
    ("SubhaloVelDisp",       (v_units, ["Subhalo", "velocity_dispersion"], None)),
    ("SubhaloVmax",          (v_units, [], None)),
    ("SubhaloVmaxRad",       (p_units, [], None)),
) + _type_fields + _sub_type_fields

class GadgetFOFFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = _particle_fields

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

    known_particle_fields = _particle_fields + \
      (("ID", ("", ["member_ids"], None)),)

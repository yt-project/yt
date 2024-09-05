from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer

m_units = "g"
p_units = "cm"
v_units = "cm / s"
r_units = "cm"

_particle_fields: KnownFieldsT = (
    ("particle_identifier", ("", [], None)),
    ("particle_position_x", (p_units, [], None)),
    ("particle_position_y", (p_units, [], None)),
    ("particle_position_z", (p_units, [], None)),
    ("particle_velocity_x", (v_units, [], None)),
    ("particle_velocity_y", (v_units, [], None)),
    ("particle_velocity_z", (v_units, [], None)),
    ("particle_mass", (m_units, [], "Virial Mass")),
    ("virial_radius", (r_units, [], "Virial Radius")),
)


class YTHaloCatalogFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = _particle_fields


class YTHaloCatalogHaloFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = _particle_fields + (("ids", ("", ["member_ids"], None)),)

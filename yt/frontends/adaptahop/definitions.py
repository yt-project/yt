"""
Data structures for AdaptaHOP



"""
import abc
from typing import Tuple, Union

from yt.funcs import mylog

ATTR_T = Tuple[Tuple[Union[Tuple[str, ...], str], int, str], ...]


def HEADER_ATTRIBUTES(*, double: bool, longint: bool) -> ATTR_T:
    int_type = "l" if longint else "i"
    float_type = "d" if double else "f"
    return (
        ("npart", 1, int_type),
        ("massp", 1, float_type),
        ("aexp", 1, float_type),
        ("omega_t", 1, float_type),
        ("age", 1, float_type),
        (("nhalos", "nsubs"), 2, "i"),
    )


ADAPTAHOP_TEMPLATES = {}


class AdaptaHOPDefTemplate(abc.ABC):
    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(*args, **kwargs)
        mylog.debug("Registering AdaptaHOP template class %s", cls.__name__)
        ADAPTAHOP_TEMPLATES[cls.__name__] = cls

    def __init__(self, longint, double_precision):
        self.longint = longint
        self.double_precision = double_precision


class AdaptaHOPOld(AdaptaHOPDefTemplate):
    @property
    def HALO_ATTRIBUTES(self) -> ATTR_T:
        int_type = "l" if self.longint else "i"
        float_type = "d" if self.double_precision else "f"
        return (
            ("npart", 1, int_type),
            ("particle_identities", -1, int_type),
            ("particle_identifier", 1, "i"),  # this is the halo id, always an int32
            ("timestep", 1, "i"),
            (
                (
                    "level",
                    "host_id",
                    "first_subhalo_id",
                    "n_subhalos",
                    "next_subhalo_id",
                ),
                5,
                "i",
            ),
            ("particle_mass", 1, float_type),
            (("raw_position_x", "raw_position_y", "raw_position_z"), 3, float_type),
            (
                ("particle_velocity_x", "particle_velocity_y", "particle_velocity_z"),
                3,
                float_type,
            ),
            (
                (
                    "particle_angular_momentum_x",
                    "particle_angular_momentum_y",
                    "particle_angular_momentum_z",
                ),
                3,
                float_type,
            ),
            (("r", "a", "b", "c"), 4, float_type),
            (("ek", "ep", "etot"), 3, float_type),
            ("spin", 1, float_type),
            (
                (
                    "virial_radius",
                    "virial_mass",
                    "virial_temperature",
                    "virial_velocity",
                ),
                4,
                float_type,
            ),
            (("rho0", "R_c"), 2, float_type),
        )


class AdaptaHOPNewNoContam(AdaptaHOPDefTemplate):
    @property
    def HALO_ATTRIBUTES(self) -> ATTR_T:
        int_type = "l" if self.longint else "i"
        float_type = "d" if self.double_precision else "f"
        return (
            ("npart", 1, int_type),
            ("particle_identities", -1, int_type),
            ("particle_identifier", 1, "i"),  # this is the halo id, always an int32
            ("timestep", 1, "i"),
            (
                (
                    "level",
                    "host_id",
                    "first_subhalo_id",
                    "n_subhalos",
                    "next_subhalo_id",
                ),
                5,
                "i",
            ),
            ("particle_mass", 1, float_type),
            ("npart_tot", 1, int_type),
            ("particle_mass_tot", 1, float_type),
            (("raw_position_x", "raw_position_y", "raw_position_z"), 3, float_type),
            (
                ("particle_velocity_x", "particle_velocity_y", "particle_velocity_z"),
                3,
                float_type,
            ),
            (
                (
                    "particle_angular_momentum_x",
                    "particle_angular_momentum_y",
                    "particle_angular_momentum_z",
                ),
                3,
                float_type,
            ),
            (("r", "a", "b", "c"), 4, float_type),
            (("ek", "ep", "etot"), 3, float_type),
            ("spin", 1, float_type),
            ("velocity_dispersion", 1, float_type),
            (
                (
                    "virial_radius",
                    "virial_mass",
                    "virial_temperature",
                    "virial_velocity",
                ),
                4,
                float_type,
            ),
            (("rmax", "vmax"), 2, float_type),
            ("concentration", 1, float_type),
            (("radius_200", "mass_200"), 2, float_type),
            (("radius_50", "mass_50"), 2, float_type),
            ("radius_profile", -1, float_type),
            ("rho_profile", -1, float_type),
            (("rho0", "R_c"), 2, float_type),
        )


class AdaptaHOPNewContam(AdaptaHOPNewNoContam):
    @property
    def HALO_ATTRIBUTES(self) -> ATTR_T:
        attrs = list(super().HALO_ATTRIBUTES)
        int_type = "l" if self.longint else "i"
        float_type = "d" if self.double_precision else "f"
        return tuple(
            attrs
            + [
                ("contaminated", 1, "i"),
                (("m_contam", "mtot_contam"), 2, float_type),
                (("n_contam", "ntot_contam"), 2, int_type),
            ]
        )

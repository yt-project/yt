"""
Data structures for AdaptaHOP



"""


HEADER_ATTRIBUTES = (
    ("npart", 1, "i"),
    ("massp", 1, "f"),
    ("aexp", 1, "f"),
    ("omega_t", 1, "f"),
    ("age", 1, "f"),
    (("nhalos", "nsubs"), 2, "i"),
)

HALO_ATTRIBUTES = (
    ("npart", 1, "i"),
    ("particle_identities", -1, "i"),
    ("particle_identifier", 1, "i"),
    ("timestep", 1, "i"),
    (("level", "host_id", "first_subhalo_id", "n_subhalos", "next_subhalo_id"), 5, "i"),
    ("particle_mass", 1, "f"),
    (("raw_position_x", "raw_position_y", "raw_position_z"), 3, "f"),
    (("particle_velocity_x", "particle_velocity_y", "particle_velocity_z"), 3, "f"),
    (
        (
            "particle_angular_momentum_x",
            "particle_angular_momentum_y",
            "particle_angular_momentum_z",
        ),
        3,
        "f",
    ),
    (("r", "a", "b", "c"), 4, "f"),
    (("ek", "ep", "etot"), 3, "f"),
    ("spin", 1, "f"),
    (("virial_radius", "virial_mass", "virial_temperature", "virial_velocity"), 4, "f"),
    (("rho0", "R_c"), 2, "f"),
)

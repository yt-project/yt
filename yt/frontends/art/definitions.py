# If not otherwise specified, we are big endian
endian = ">"

fluid_fields = [
    "Density",
    "TotalEnergy",
    "XMomentumDensity",
    "YMomentumDensity",
    "ZMomentumDensity",
    "Pressure",
    "Gamma",
    "GasEnergy",
    "MetalDensitySNII",
    "MetalDensitySNIa",
    "PotentialNew",
    "PotentialOld",
]

hydro_struct = [("pad1", ">i"), ("idc", ">i"), ("iOctCh", ">i")]
for field in fluid_fields:
    hydro_struct += ((field, ">f"),)
hydro_struct += (("pad2", ">i"),)

particle_fields = [
    "particle_mass",  # stars have variable mass
    "particle_index",
    "particle_type",
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_velocity_x",
    "particle_velocity_y",
    "particle_velocity_z",
    "particle_mass_initial",
    "particle_creation_time",
    "particle_metallicity1",
    "particle_metallicity2",
    "particle_metallicity",
]

particle_star_fields = [
    "particle_mass",
    "particle_mass_initial",
    "particle_creation_time",
    "particle_metallicity1",
    "particle_metallicity2",
    "particle_metallicity",
]


filename_pattern = {
    "amr": ["10MpcBox_", ".d"],
    "particle_header": ["PMcrd", ".DAT"],
    "particle_data": ["PMcrs", ".DAT"],
    "particle_stars": ["stars", ".dat"],
}

amr_header_struct = [
    ("jname", 1, "256s"),
    (("istep", "t", "dt", "aexpn", "ainit"), 1, "iddff"),
    (("boxh", "Om0", "Oml0", "Omb0", "hubble"), 5, "f"),
    ("nextras", 1, "i"),
    (("extra1", "extra2"), 2, "f"),
    ("lextra", 1, "512s"),
    (("min_level", "max_level"), 2, "i"),
]

particle_header_struct = [
    (
        (
            "header",
            "aexpn",
            "aexp0",
            "amplt",
            "astep",
            "istep",
            "partw",
            "tintg",
            "Ekin",
            "Ekin1",
            "Ekin2",
            "au0",
            "aeu0",
            "Nrow",
            "Ngridc",
            "Nspecies",
            "Nseed",
            "Om0",
            "Oml0",
            "hubble",
            "Wp5",
            "Ocurv",
            "Omb0",
            "extras",
            "unknown",
        ),
        1,
        "45sffffi" + "fffffff" + "iiii" + "ffffff" + "396s" + "f",
    )
]

dmparticle_header_struct = [
    (
        "header",
        "aexpn",
        "aexp0",
        "amplt",
        "astep",
        "istep",
        "partw",
        "tintg",
        "Ekin",
        "Ekin1",
        "Ekin2",
        "au0",
        "aeu0",
        "Nrow",
        "Ngridc",
        "Nspecies",
        "Nseed",
        "Om0",
        "Oml0",
        "hubble",
        "Wp5",
        "Ocurv",
        "wspecies",
        "lspecies",
        "extras",
        "boxsize",
    ),
    (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 79, 1),
]

star_struct = [
    (">d", ("t_stars", "a_stars")),
    (">i", "nstars"),
    (">d", ("ws_old", "ws_oldi")),
    (">f", "particle_mass"),
    (">f", "particle_mass_initial"),
    (">f", "particle_creation_time"),
    (">f", "particle_metallicity1"),
    (">f", "particle_metallicity2"),
]

star_name_map = {
    "particle_mass": "mass",
    "particle_mass_initial": "imass",
    "particle_creation_time": "tbirth",
    "particle_metallicity1": "metallicity1",
    "particle_metallicity2": "metallicity2",
    "particle_metallicity": "metallicity",
}

constants = {
    "Y_p": 0.245,
    "gamma": 5.0 / 3.0,
    "T_CMB0": 2.726,
    "T_min": 300.0,
    "wmu": 4.0 / (8.0 - 5.0 * 0.245),
}

seek_extras = 137

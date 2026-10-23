import re

# RAMSES II (and mini-ramses) outputs live in a folder named output_XXXXX
# (or backup_XXXXX for restart dumps) containing
# - info.txt, params.bin, namelist.txt, compilation.txt, timer.txt
# - amr.NNNNN                          the oct structure
# - hydro.NNNNN, grav.NNNNN, rt.NNNNN, cr.NNNNN
#                                      cell-based fields (one block per oct)
# - part.NNNNN, star.NNNNN, sink.NNNNN, tree.NNNNN, trac.NNNNN, dust.NNNNN,
#   core.NNNNN                         one file per particle type
# - <prefix>_header.txt                text descriptors of the binary files
# where NNNNN is the file index, running from 1 to nfile.
# The layouts follow utils/py/ramses_io.py in the RAMSES II repository.
OUTPUT_DIR_RE = re.compile(r"output_(\d{5})")

FLUID_PREFIXES = ("hydro", "grav", "rt", "cr")
PARTICLE_PREFIXES = ("part", "star", "sink", "tree", "trac", "dust", "core")

RAMSES2_FILE_RE = re.compile(
    rf"(({'|'.join(('amr',) + FLUID_PREFIXES + PARTICLE_PREFIXES)})\.\d{{5}}"
    rf"|({'|'.join(FLUID_PREFIXES + PARTICLE_PREFIXES)})_header\.txt"
    r"|info\.txt|params\.bin)"
)

# yt field type of each fluid file. Fields keep the names of
# <prefix>_header.txt (density, velocity_x, thermal_pressure, ...).
FLUID_TYPES = {
    "hydro": "ramses",
    "grav": "gravity",
    "rt": "rt",
    "cr": "cr",
}

# yt name of the particle fields, without the _x, _y, _z suffix of vector
# fields. Other fields are named "particle_" + their name in the header.
PARTICLE_FIELD_NAMES = {
    "pos": "particle_position",
    "vel": "particle_velocity",
    "accel": "particle_acceleration",
    "angmom": "particle_angular_momentum",
    "mass": "particle_mass",
    "birth_id": "particle_identity",
    "level": "particle_level",
    "birth_date": "particle_birth_time",
    "merging_date": "particle_merging_time",
    "metallicity": "particle_metallicity",
}

# Content of particle files written before the particle headers listed the
# type and number of components of each field (mini-ramses and RAMSES II
# before October 2026). "id" fields are 4 or 8 bytes long depending on the
# compilation flags of the code: their size is inferred from the file size.
LEGACY_PARTICLE_FIELDS = {
    "pos": ("f4", 3),
    "vel": ("f4", 3),
    "mass": ("f4", 1),
    "potential": ("f4", 1),
    "metallicity": ("f4", 1),
    "accel": ("f4", 3),
    "angmom": ("f4", 3),
    "birth_date": ("f4", 1),
    "merging_date": ("f4", 1),
    "level": ("i4", 1),
    "birth_id": ("id", 1),
    "merging_id": ("id", 1),
    "tracking_id": ("id", 1),
}

field_aliases = {
    "standard_five": (
        "density",
        "velocity_x",
        "velocity_y",
        "velocity_z",
        "thermal_pressure",
    ),
}

import re

from yt.funcs import mylog

# mini-ramses uses different file naming than RAMSES:
# - output_XXXXX/info.txt (no number in info filename)
# - output_XXXXX/amr.NNNNN (dot-separated, no "out" suffix)
# - output_XXXXX/hydro.NNNNN
# - output_XXXXX/part.NNNNN, star.NNNNN, sink.NNNNN, tree.NNNNN
# - output_XXXXX/grav.NNNNN
# - output_XXXXX/params.bin
# - output_XXXXX/hydro_header.txt, part_header.txt, etc.

OUTPUT_DIR_EXP = r"output_(\d{5})"
OUTPUT_DIR_RE = re.compile(OUTPUT_DIR_EXP)

# Matches mini-ramses file patterns (dot-separated, no "out" suffix)
# e.g. amr.00001, hydro.00001, part.00001, info.txt
MINI_RAMSES_FILE_RE = re.compile(
    r"((amr|hydro|part|grav|star|sink|tree|trac|rt)\.\d{5}"
    r"|info\.txt|params\.bin)"
)

# Regular expressions used to parse file descriptors
VERSION_RE = re.compile(r"# version: *(\d+)")
VAR_DESC_RE = re.compile(r"\s*([^\s]+),\s*([^\s]+),\s*([^\s]+)")

field_aliases = {
    "standard_five": (
        "Density",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Pressure",
    ),
    "standard_six": (
        "Density",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Pressure",
        "Metallicity",
    ),
}

# mini-ramses uses separate files for each particle type,
# so the family concept is simpler
particle_families = {
    "DM": 1,
    "star": 2,
    "cloud": 3,
}

# These functions are RAMSES-specific
import re

from yt.funcs import mylog
from yt.utilities.configure import YTConfig, configuration_callbacks


def ramses_header(hvals):
    header = (
        ("ncpu", 1, "i"),
        ("ndim", 1, "i"),
        ("nx", 3, "i"),
        ("nlevelmax", 1, "i"),
        ("ngridmax", 1, "i"),
        ("nboundary", 1, "i"),
        ("ngrid_current", 1, "i"),
        ("boxlen", 1, "d"),
        ("nout", 3, "i"),
    )
    yield header
    # TODO: REMOVE
    noutput, iout, ifout = hvals["nout"]
    next_set = (
        ("tout", noutput, "d"),
        ("aout", noutput, "d"),
        ("t", 1, "d"),
        ("dtold", hvals["nlevelmax"], "d"),
        ("dtnew", hvals["nlevelmax"], "d"),
        ("nstep", 2, "i"),
        ("stat", 3, "d"),
        ("cosm", 7, "d"),
        ("timing", 5, "d"),
        ("mass_sph", 1, "d", True),
    )
    yield next_set


field_aliases = {
    "standard_five": ("Density", "x-velocity", "y-velocity", "z-velocity", "Pressure"),
    "standard_six": (
        "Density",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Pressure",
        "Metallicity",
    ),
}

## Regular expressions used to parse file descriptors
VERSION_RE = re.compile(r"# version: *(\d+)")
# This will match comma-separated strings, discarding whitespaces
# on the left hand side
VAR_DESC_RE = re.compile(r"\s*([^\s]+),\s*([^\s]+),\s*([^\s]+)")

OUTPUT_DIR_EXP = r"output_(\d{5})"
OUTPUT_DIR_RE = re.compile(OUTPUT_DIR_EXP)
STANDARD_FILE_RE = re.compile(r"((amr|hydro|part|grav)_\d{5}\.out\d{5}|info_\d{5}.txt)")


## Configure family mapping
particle_families = {
    "DM": 1,
    "star": 2,
    "cloud": 3,
    "dust": 4,
    "star_tracer": -2,
    "cloud_tracer": -3,
    "dust_tracer": -4,
    "gas_tracer": 0,
}


def _setup_ramses_particle_families(ytcfg: YTConfig) -> None:
    if not ytcfg.has_section("ramses-families"):
        return
    for key in particle_families.keys():
        val = ytcfg.get("ramses-families", key, callback=None)
        if val is not None:
            mylog.info(
                "Changing family %s from %s to %s", key, particle_families[key], val
            )
            particle_families[key] = val


configuration_callbacks.append(_setup_ramses_particle_families)

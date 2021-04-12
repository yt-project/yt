from .physical_ratios import (
    au_per_mpc,
    cm_per_mpc,
    km_per_mpc,
    kpc_per_mpc,
    miles_per_mpc,
    mpc_per_mpc,
    pc_per_mpc,
    rsun_per_mpc,
    sec_per_day,
    sec_per_Gyr,
    sec_per_Myr,
    sec_per_year,
)

# The number of levels we expect to have at most
MAXLEVEL = 48

# How many of each thing are in an Mpc
mpc_conversion = {
    "Mpc": mpc_per_mpc,
    "mpc": mpc_per_mpc,
    "kpc": kpc_per_mpc,
    "pc": pc_per_mpc,
    "au": au_per_mpc,
    "rsun": rsun_per_mpc,
    "miles": miles_per_mpc,
    "km": km_per_mpc,
    "cm": cm_per_mpc,
}

# Nicely formatted versions of common length units
formatted_length_unit_names = {
    "au": "AU",
    "rsun": r"R_\odot",
    "code_length": r"code\ length",
}

# How many seconds are in each thing
sec_conversion = {
    "Gyr": sec_per_Gyr,
    "Myr": sec_per_Myr,
    "years": sec_per_year,
    "days": sec_per_day,
}

from yt._maintenance.deprecation import issue_deprecation_warning, warnings


def boxlib_deprecation():
    warnings.simplefilter("always")
    issue_deprecation_warning(
        "The historic 'boxlib' frontend is \n"
        "deprecated as it has been renamed 'amrex'. "
        "Existing and future work should instead reference the 'amrex' frontend.",
        stacklevel=3,
        since="TBD",
    )
    warnings.resetwarnings()

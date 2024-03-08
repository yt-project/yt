from ...amrex import tests
from yt._maintenance.deprecation import issue_deprecation_warning

issue_deprecation_warning(
    "The historic 'boxlib' frontend is \n"
    "deprecated as it has been renamed 'amrex'. "
    "Future work should reference the 'amrex' frontend.",
    stacklevel=3,
    since="TBD",
)

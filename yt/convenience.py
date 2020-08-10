# this is a deprecated module
from .funcs import issue_deprecation_warning
from .loaders import load, load_simulation, simulation

issue_deprecation_warning(
    "importing from yt.convenience is deprecated in favor of yt.loaders\n"
    "This will become an error in yt 4.1"
)
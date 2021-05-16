import sys
from importlib import import_module
from typing import Optional

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


class NotAModule:
    """
    A class to implement an informative error message that will be outputted if
    someone tries to use an on-demand import without having the requisite
    package installed.
    """

    def __init__(self, pkg_name):
        self.pkg_name = pkg_name
        self.error = ImportError(
            f"This functionality requires the {self.pkg_name} package to be installed."
        )

    def __getattr__(self, item):
        raise self.error

    def __call__(self, *args, **kwargs):
        raise self.error


class OnDemandImport:
    """
    Implement a safely importable drop-in replacement for a module that
    might not be installed, and fail *only* on the first attempt to use it.
    """

    def __init__(self, name: str, instructions: Optional[str] = None):
        self._name = name
        self._submodules = {}
        msg = f"Failed to import an optional dependency from yt ({self._name})."
        if instructions:
            msg += f"\n{instructions}"
        self._error = ImportError(msg)

    @cached_property
    def module(self):
        try:
            return import_module(self._name)
        except ImportError as exc:
            raise self._error from exc

    def __getattr__(self, key):
        # the recursive nature of this method is intended to support
        # arbitrary levels of nesting for submodules that need to be
        # explicitly imported on their own, like scipy.spatial or astropy.io.fits
        if key in self._submodules:
            return self._submodules[key]
        try:
            return getattr(self.module, key)
        except AttributeError:
            self._submodules[key] = OnDemandImport(f"{self._name}.{key}")
            return self._submodules[key]


if any(s in sys.version for s in ("Anaconda", "Continuum")):
    # conda distributions of cartopy don't have issues with the GEOS library
    # so we can keep the error message short
    # see https://github.com/yt-project/yt/pull/1966
    _cartopy_instructions = None
else:
    _cartopy_instructions = (
        "Try installing proj4 and geos with your package manager "
        "and building shapely and cartopy from source via Pypi with:\n"
        "`python -m pip install --no-binary :all: shapely cartopy`\n"
        "For further instruction please refer to the yt documentation."
    )

optional_dependencies = {
    "astropy": None,
    "bottle": None,
    "cartopy": _cartopy_instructions,
    "Firefly": None,
    "f90nml": None,
    # Breaking alphabetical ordering on purpose here.
    # Importing h5py before netCDF4 can lead to file lock errors
    # on some systems (including travis-ci)
    # see https://github.com/pydata/xarray/issues/2560
    "netCDF4": None,
    "h5py": None,
    "libconf": None,
    "miniball": "This package can be installed from PyPI with\n`python3 -m pip install MiniballCpp`",
    "mpi4py": None,
    "nose": None,
    "pandas": None,
    "pooch": None,
    "requests": None,
    "scipy": None,
    "yaml": None,
}

_module = sys.modules[__name__]
for package_name, instructions in optional_dependencies.items():
    setattr(_module, package_name, OnDemandImport(package_name, instructions))

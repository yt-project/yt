import sys
from importlib import import_module
from typing import Optional


class NotAModule:
    """
    A class to implement an informative error message that will be outputted if
    someone tries to use an on-demand import without having the requisite
    module installed.
    """

    def __init__(self, pkg_name):
        self.pkg_name = pkg_name
        self.error = ImportError(
            f"This functionality requires the {self.pkg_name} module to be installed."
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

    def __init__(self, name: str, message: Optional[str] = None):
        self._name = name
        self._module = None
        self._submodules = {}
        msg = f"This functionality requires the {self._name} module to be installed."
        if message is not None:
            msg += f" {message}"
        self._error = ImportError(msg)

    @property
    def module(self):
        if self._module is None:
            try:
                self._module = import_module(self._name)
            except ImportError as exc:
                raise self._error from exc
        return self._module

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


targets = [
    "astropy",
    "bottle",
    # "f90nml",
    "firefly_api",
    "h5py",
    "libconf",
    "mpi4py",
    "netCDF4",
    "nose",
    "pandas",
    "pooch",
    "requests",
    "resource",
    "scipy",
    "yaml",
    "yt_astro_analysis",
]
_module = sys.modules[__name__]
for target in targets:
    setattr(_module, target, OnDemandImport(target))

_cartopy_message = ""
if not any(s in sys.version for s in ("Anaconda", "Continuum")):
    # the conda-based installs of cartopy don't have issues with the
    # GEOS library, so the error message for users with conda can be
    # relatively short. Discussion related to this is in
    # yt-project/yt#1966
    _cartopy_message = (
        "Try installing proj4 and "
        "geos with your package manager and building shapely "
        "and cartopy from source with: \n\n "
        "pip install --no-binary :all: shapely cartopy \n\n"
        "For further instruction please refer to the "
        "yt documentation."
    )

cartopy = OnDemandImport("cartopy", _cartopy_message)

miniball = OnDemandImport(
    "miniball",
    (
        "Installation instructions can be found at "
        "https://github.com/weddige/miniball or alternatively you can "
        "install via `pip install MiniballCpp`."
    ),
)

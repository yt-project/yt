import sys
from functools import wraps
from importlib.util import find_spec
from typing import Type


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

    def __repr__(self) -> str:
        return f"NotAModule({self.pkg_name!r})"


class OnDemand:
    _default_factory: Type[NotAModule] = NotAModule

    def __init_subclass__(cls):
        if not cls.__name__.endswith("_imports"):
            raise TypeError(f"class {cls}'s name needs to be suffixed '_imports'")

    def __new__(cls):
        if cls is OnDemand:
            raise TypeError("The OnDemand base class cannot be instanciated.")
        else:
            return object.__new__(cls)

    @property
    def _name(self) -> str:
        _name, _, _suffix = self.__class__.__name__.rpartition("_")
        return _name

    @property
    def __is_available__(self) -> bool:
        # special protocol to support testing framework
        return find_spec(self._name) is not None


def safe_import(func):
    @property
    @wraps(func)
    def inner(self):
        try:
            return func(self)
        except ImportError:
            return self._default_factory(self._name)

    return inner


class netCDF4_imports(OnDemand):
    def __init__(self):
        # this ensures the import ordering between netcdf4 and h5py. If h5py is
        # imported first, can get file lock errors on some systems (including travis-ci)
        # so we need to do this before initializing h5py_imports()!
        # similar to this issue https://github.com/pydata/xarray/issues/2560
        try:
            import netCDF4  # noqa F401
        except ImportError:
            pass

    @safe_import
    def Dataset(self):
        from netCDF4 import Dataset

        return Dataset


_netCDF4 = netCDF4_imports()


class astropy_imports(OnDemand):
    @safe_import
    def log(self):
        from astropy import log

        if log.exception_logging_enabled():
            log.disable_exception_logging()

        return log

    @safe_import
    def pyfits(self):
        from astropy.io import fits

        return fits

    @safe_import
    def pywcs(self):
        import astropy.wcs as pywcs

        self.log
        return pywcs

    @safe_import
    def units(self):
        from astropy import units

        self.log
        return units

    @safe_import
    def conv(self):
        import astropy.convolution as conv

        self.log
        return conv

    @safe_import
    def time(self):
        import astropy.time as time

        self.log
        return time

    @safe_import
    def wcsaxes(self):
        from astropy.visualization import wcsaxes

        self.log
        return wcsaxes


_astropy = astropy_imports()


class NotCartopy(NotAModule):
    """
    A custom class to return error messages dependent on system installation
    for cartopy imports.
    """

    def __init__(self, pkg_name):
        self.pkg_name = pkg_name
        if any(s in sys.version for s in ("Anaconda", "Continuum")):
            # the conda-based installs of cartopy don't have issues with the
            # GEOS library, so the error message for users with conda can be
            # relatively short. Discussion related to this is in
            # yt-project/yt#1966
            self.error = ImportError(
                "This functionality requires the %s "
                "package to be installed." % self.pkg_name
            )
        else:
            self.error = ImportError(
                f"This functionality requires the {pkg_name} "
                "package to be installed.\n"
                "For further instruction please refer to Cartopy's documentation\n"
                "https://scitools.org.uk/cartopy/docs/latest/installing.html"
            )


class cartopy_imports(OnDemand):

    _default_factory = NotCartopy

    @safe_import
    def crs(self):
        import cartopy.crs as crs

        return crs


_cartopy = cartopy_imports()


class pooch_imports(OnDemand):
    @safe_import
    def HTTPDownloader(self):
        from pooch import HTTPDownloader

        return HTTPDownloader

    @safe_import
    def utils(self):
        from pooch import utils

        return utils

    @safe_import
    def create(self):
        from pooch import create

        return create


_pooch = pooch_imports()


class pyart_imports(OnDemand):
    @safe_import
    def io(self):
        from pyart import io

        return io

    @safe_import
    def map(self):
        from pyart import map

        return map


_pyart = pyart_imports()


class xarray_imports(OnDemand):
    @safe_import
    def open_dataset(self):
        from xarray import open_dataset

        return open_dataset


_xarray = xarray_imports()


class scipy_imports(OnDemand):
    @safe_import
    def signal(self):
        from scipy import signal

        return signal

    @safe_import
    def spatial(self):
        from scipy import spatial

        return spatial

    @safe_import
    def ndimage(self):
        from scipy import ndimage

        return ndimage

    # Optimize is not presently used by yt, but appears here to enable
    # required functionality in yt extension, trident

    @safe_import
    def optimize(self):
        from scipy import optimize

        return optimize


_scipy = scipy_imports()


class h5py_imports(OnDemand):
    @safe_import
    def File(self):
        from h5py import File

        return File

    @safe_import
    def Group(self):
        from h5py import Group

        return Group

    @safe_import
    def Dataset(self):
        from h5py import Dataset

        return Dataset

    @safe_import
    def get_config(self):
        from h5py import get_config

        return get_config

    @safe_import
    def h5f(self):
        from h5py import h5f

        return h5f

    @safe_import
    def h5p(self):
        from h5py import h5p

        return h5p

    @safe_import
    def h5d(self):
        from h5py import h5d

        return h5d

    @safe_import
    def h5s(self):
        from h5py import h5s

        return h5s


_h5py = h5py_imports()


class nose_imports(OnDemand):
    @safe_import
    def run(self):
        from nose import run

        return run


_nose = nose_imports()


class libconf_imports(OnDemand):
    @safe_import
    def load(self):
        from libconf import load

        return load


_libconf = libconf_imports()


class yaml_imports(OnDemand):
    @safe_import
    def load(self):
        from yaml import load

        return load

    @safe_import
    def FullLoader(self):
        from yaml import FullLoader

        return FullLoader


_yaml = yaml_imports()


class NotMiniball(NotAModule):
    def __init__(self, pkg_name):
        super().__init__(pkg_name)
        str = (
            "This functionality requires the %s package to be installed. "
            "Installation instructions can be found at "
            "https://github.com/weddige/miniball or alternatively you can "
            "install via `python -m pip install MiniballCpp`."
        )
        self.error = ImportError(str % self.pkg_name)


class miniball_imports(OnDemand):
    @safe_import
    def Miniball(self):
        from miniball import Miniball

        return Miniball


_miniball = miniball_imports()


class f90nml_imports(OnDemand):
    @safe_import
    def read(self):
        from f90nml import read

        return read

    @safe_import
    def Namelist(self):
        from f90nml import Namelist

        return Namelist


_f90nml = f90nml_imports()


class requests_imports(OnDemand):
    @safe_import
    def post(self):
        from requests import post

        return post

    @safe_import
    def put(self):
        from requests import put

        return put

    @safe_import
    def codes(self):
        from requests import codes

        return codes

    @safe_import
    def get(self):
        from requests import get

        return get

    @safe_import
    def exceptions(self):
        from requests import exceptions

        return exceptions


_requests = requests_imports()


class pandas_imports(OnDemand):
    @safe_import
    def NA(self):
        from pandas import NA

        return NA

    @safe_import
    def DataFrame(self):
        from pandas import DataFrame

        return DataFrame

    @safe_import
    def concat(self):
        from pandas import concat

        return concat


_pandas = pandas_imports()


class firefly_imports(OnDemand):
    @safe_import
    def data_reader(self):
        import firefly.data_reader as data_reader

        return data_reader

    @safe_import
    def server(self):
        import firefly.server as server

        return server


_firefly = firefly_imports()


# Note: ratarmount may fail with an OSError on import if libfuse is missing
class ratarmount_imports(OnDemand):
    @safe_import
    def TarMount(self):
        from ratarmount import TarMount

        return TarMount

    @safe_import
    def fuse(self):
        from ratarmount import fuse

        return fuse


_ratarmount = ratarmount_imports()

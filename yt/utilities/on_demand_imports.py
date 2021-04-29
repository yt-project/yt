import sys
from distutils.version import LooseVersion


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
                "This functionality requires the %s "
                "package to be installed. Try installing proj4 and "
                "geos with your package manager and building shapely "
                "and cartopy from source with: \n \n "
                "pip install --no-binary :all: shapely cartopy \n \n"
                "For further instruction please refer to the "
                "yt documentation." % self.pkg_name
            )


class netCDF4_imports:
    _name = "netCDF4"
    _Dataset = None

    def __init__(self):
        # this ensures the import ordering between netcdf4 and h5py. If h5py is
        # imported first, can get file lock errors on some systems (including travis-ci)
        # so we need to do this before initializing h5py_imports()!
        # similar to this issue https://github.com/pydata/xarray/issues/2560
        try:
            import netCDF4  # noqa F401
        except ImportError:
            pass

    @property
    def Dataset(self):
        if self._Dataset is None:
            try:
                from netCDF4 import Dataset
            except ImportError:
                Dataset = NotAModule(self._name)
            self._Dataset = Dataset
        return self._Dataset


_netCDF4 = netCDF4_imports()


class astropy_imports:
    _name = "astropy"
    _pyfits = None

    @property
    def pyfits(self):
        if self._pyfits is None:
            try:
                import astropy.io.fits as pyfits

                self.log
            except ImportError:
                pyfits = NotAModule(self._name)
            self._pyfits = pyfits
        return self._pyfits

    _pywcs = None

    @property
    def pywcs(self):
        if self._pywcs is None:
            try:
                import astropy.wcs as pywcs

                self.log
            except ImportError:
                pywcs = NotAModule(self._name)
            self._pywcs = pywcs
        return self._pywcs

    _log = None

    @property
    def log(self):
        if self._log is None:
            try:
                from astropy import log

                if log.exception_logging_enabled():
                    log.disable_exception_logging()
            except ImportError:
                log = NotAModule(self._name)
            self._log = log
        return self._log

    _units = None

    @property
    def units(self):
        if self._units is None:
            try:
                from astropy import units

                self.log
            except ImportError:
                units = NotAModule(self._name)
            self._units = units
        return self._units

    _conv = None

    @property
    def conv(self):
        if self._conv is None:
            try:
                import astropy.convolution as conv

                self.log
            except ImportError:
                conv = NotAModule(self._name)
            self._conv = conv
        return self._conv

    _time = None

    @property
    def time(self):
        if self._time is None:
            try:
                import astropy.time as time

                self.log
            except ImportError:
                time = NotAModule(self._name)
            self._time = time
        return self._time

    _wcsaxes = None

    @property
    def wcsaxes(self):
        if self._wcsaxes is None:
            try:
                import astropy.visualization.wcsaxes as wcsaxes

                self.log
            except ImportError:
                wcsaxes = NotAModule(self._name)
            self._wcsaxes = wcsaxes
        return self._wcsaxes

    _version = None

    @property
    def __version__(self):
        if self._version is None:
            try:
                import astropy

                version = astropy.__version__
            except ImportError:
                version = NotAModule(self._name)
            self._version = version
        return self._version


_astropy = astropy_imports()


class cartopy_imports:
    _name = "cartopy"

    _crs = None

    @property
    def crs(self):
        if self._crs is None:
            try:
                import cartopy.crs as crs
            except ImportError:
                crs = NotCartopy(self._name)
            self._crs = crs
        return self._crs

    _version = None

    @property
    def __version__(self):
        if self._version is None:
            try:
                import cartopy

                version = cartopy.__version__
            except ImportError:
                version = NotCartopy(self._name)
            self._version = version
        return self._version


_cartopy = cartopy_imports()


class pooch_imports:
    _name = "pooch"
    _module = None

    def __init__(self):
        try:
            import pooch as myself

            self._module = myself
        except ImportError:
            self._module = NotAModule(self._name)

    def __getattr__(self, attr):
        return getattr(self._module, attr)


_pooch = pooch_imports()


class scipy_imports:
    _name = "scipy"
    _integrate = None

    @property
    def integrate(self):
        if self._integrate is None:
            try:
                import scipy.integrate as integrate
            except ImportError:
                integrate = NotAModule(self._name)
            self._integrate = integrate
        return self._integrate

    _stats = None

    @property
    def stats(self):
        if self._stats is None:
            try:
                import scipy.stats as stats
            except ImportError:
                stats = NotAModule(self._name)
            self._stats = stats
        return self._stats

    _optimize = None

    @property
    def optimize(self):
        if self._optimize is None:
            try:
                import scipy.optimize as optimize
            except ImportError:
                optimize = NotAModule(self._name)
            self._optimize = optimize
        return self._optimize

    _interpolate = None

    @property
    def interpolate(self):
        if self._interpolate is None:
            try:
                import scipy.interpolate as interpolate
            except ImportError:
                interpolate = NotAModule(self._name)
            self._interpolate = interpolate
        return self._interpolate

    _special = None

    @property
    def special(self):
        if self._special is None:
            try:
                import scipy.special as special
            except ImportError:
                special = NotAModule(self._name)
            self._special = special
        return self._special

    _signal = None

    @property
    def signal(self):
        if self._signal is None:
            try:
                import scipy.signal as signal
            except ImportError:
                signal = NotAModule(self._name)
            self._signal = signal
        return self._signal

    _spatial = None

    @property
    def spatial(self):
        if self._spatial is None:
            try:
                import scipy.spatial as spatial
            except ImportError:
                spatial = NotAModule(self._name)
            self._spatial = spatial
        return self._spatial

    _ndimage = None

    @property
    def ndimage(self):
        if self._ndimage is None:
            try:
                import scipy.ndimage as ndimage
            except ImportError:
                ndimage = NotAModule(self._name)
            self._ndimage = ndimage
        return self._ndimage


_scipy = scipy_imports()


class h5py_imports:
    _name = "h5py"
    _err = None

    def __init__(self):
        try:
            import h5py

            if LooseVersion(h5py.__version__) < LooseVersion("2.4.0"):
                self._err = RuntimeError(
                    "yt requires h5py version 2.4.0 or newer, "
                    'please update h5py with e.g. "pip install -U h5py" '
                    "and try again"
                )
        except ImportError:
            pass
        super().__init__()

    _File = None

    @property
    def File(self):
        if self._err:
            raise self._err
        if self._File is None:
            try:
                from h5py import File
            except ImportError:
                File = NotAModule(self._name)
            self._File = File
        return self._File

    _Group = None

    @property
    def Group(self):
        if self._err:
            raise self._err
        if self._Group is None:
            try:
                from h5py import Group
            except ImportError:
                Group = NotAModule(self._name)
            self._Group = Group
        return self._Group

    _Dataset = None

    @property
    def Dataset(self):
        if self._err:
            raise self._err
        if self._Dataset is None:
            try:
                from h5py import Dataset
            except ImportError:
                Dataset = NotAModule(self._name)
            self._Dataset = Dataset
        return self._Dataset

    ___version__ = None

    @property
    def __version__(self):
        if self._err:
            raise self._err
        if self.___version__ is None:
            try:
                from h5py import __version__
            except ImportError:
                __version__ = NotAModule(self._name)
            self.___version__ = __version__
        return self.___version__

    _get_config = None

    @property
    def get_config(self):
        if self._err:
            raise self._err
        if self._get_config is None:
            try:
                from h5py import get_config
            except ImportError:
                get_config = NotAModule(self._name)
            self._get_config = get_config
        return self._get_config

    _h5f = None

    @property
    def h5f(self):
        if self._err:
            raise self._err
        if self._h5f is None:
            try:
                import h5py.h5f as h5f
            except ImportError:
                h5f = NotAModule(self._name)
            self._h5f = h5f
        return self._h5f

    _h5p = None

    @property
    def h5p(self):
        if self._err:
            raise self._err
        if self._h5p is None:
            try:
                import h5py.h5p as h5p
            except ImportError:
                h5p = NotAModule(self._name)
            self._h5p = h5p
        return self._h5p

    _h5d = None

    @property
    def h5d(self):
        if self._err:
            raise self._err
        if self._h5d is None:
            try:
                import h5py.h5d as h5d
            except ImportError:
                h5d = NotAModule(self._name)
            self._h5d = h5d
        return self._h5d

    _h5s = None

    @property
    def h5s(self):
        if self._err:
            raise self._err
        if self._h5s is None:
            try:
                import h5py.h5s as h5s
            except ImportError:
                h5s = NotAModule(self._name)
            self._h5s = h5s
        return self._h5s

    _version = None

    @property
    def version(self):
        if self._err:
            raise self._err
        if self._version is None:
            try:
                import h5py.version as version
            except ImportError:
                version = NotAModule(self._name)
            self._version = version
        return self._version


_h5py = h5py_imports()


class nose_imports:
    _name = "nose"
    _run = None

    @property
    def run(self):
        if self._run is None:
            try:
                from nose import run
            except ImportError:
                run = NotAModule(self._name)
            self._run = run
        return self._run


_nose = nose_imports()


class libconf_imports:
    _name = "libconf"
    _load = None

    @property
    def load(self):
        if self._load is None:
            try:
                from libconf import load
            except ImportError:
                load = NotAModule(self._name)
            self._load = load
        return self._load


_libconf = libconf_imports()


class yaml_imports:
    _name = "yaml"
    _load = None
    _FullLoader = None

    @property
    def load(self):
        if self._load is None:
            try:
                from yaml import load
            except ImportError:
                load = NotAModule(self._name)
            self._load = load
        return self._load

    @property
    def FullLoader(self):
        if self._FullLoader is None:
            try:
                from yaml import FullLoader
            except ImportError:
                FullLoader = NotAModule(self._name)
            self._FullLoader = FullLoader
        return self._FullLoader


_yaml = yaml_imports()


class NotMiniball(NotAModule):
    def __init__(self, pkg_name):
        super().__init__(pkg_name)
        str = (
            "This functionality requires the %s package to be installed. "
            "Installation instructions can be found at "
            "https://github.com/weddige/miniball or alternatively you can "
            "install via `pip install MiniballCpp`."
        )
        self.error = ImportError(str % self.pkg_name)


class miniball_imports:
    _name = "miniball"
    _Miniball = None

    @property
    def Miniball(self):
        if self._Miniball is None:
            try:
                from miniball import Miniball
            except ImportError:
                Miniball = NotMiniball(self._name)
            self._Miniball = Miniball
        return self._Miniball


_miniball = miniball_imports()


class f90nml_imports:
    _name = "f90nml"
    _module = None

    def __init__(self):
        try:
            import f90nml as myself

            self._module = myself
        except ImportError:
            self._module = NotAModule(self._name)

    def __getattr__(self, attr):
        return getattr(self._module, attr)


_f90nml = f90nml_imports()


class requests_imports:
    _name = "requests"
    _module = None

    def __init__(self):
        try:
            import requests as myself

            self._module = myself
        except ImportError:
            self._module = NotAModule(self._name)

    def __getattr__(self, attr):
        return getattr(self._module, attr)


_requests = requests_imports()


class pandas_imports:
    _name = "pandas"
    _module = None

    def __init__(self):
        try:
            import pandas as myself

            self._module = myself
        except ImportError:
            self._module = NotAModule(self._name)

    def __getattr__(self, attr):
        return getattr(self._module, attr)


_pandas = pandas_imports()

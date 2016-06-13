"""
A set of convenient on-demand imports
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from pkg_resources import parse_version

class NotAModule(object):
    """
    A class to implement an informative error message that will be outputted if
    someone tries to use an on-demand import without having the requisite
    package installed.
    """
    def __init__(self, pkg_name):
        self.pkg_name = pkg_name
        self.error = ImportError(
            "This functionality requires the %s "
            "package to be installed." % self.pkg_name)

    def __getattr__(self, item):
        raise self.error

    def __call__(self, *args, **kwargs):
        raise self.error

class astropy_imports(object):
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

_astropy = astropy_imports()

class scipy_imports(object):
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

_scipy = scipy_imports()

class h5py_imports(object):
    _name = "h5py"
    _err = None

    def __init__(self):
        try:
            import h5py
            if parse_version(h5py.__version__) < parse_version('2.4.0'):
                self._err = RuntimeError(
                    'yt requires h5py version 2.4.0 or newer, '
                    'please update h5py with e.g. "pip install -U h5py" '
                    'and try again')
        except ImportError:
            pass
        super(h5py_imports, self).__init__()

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

class nose_imports(object):
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

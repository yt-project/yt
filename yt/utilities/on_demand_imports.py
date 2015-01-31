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

class NotAModule(object):
    """
    A class to implement an informative error message that will be outputted if
    someone tries to use an on-demand import without having the requisite package installed.
    """
    def __init__(self, pkg_name):
        self.pkg_name = pkg_name

    def __getattr__(self, item):
        raise ImportError("This functionality requires the %s package to be installed."
                          % self.pkg_name)

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

_astropy = astropy_imports()

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

_scipy = scipy_imports()
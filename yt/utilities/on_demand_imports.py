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
    _pyfits = None
    @property
    def pyfits(self):
        if self._pyfits is None:
            try:
                import astropy.io.fits as pyfits
                self.log
            except ImportError:
                pyfits = NotAModule("astropy")
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
                pywcs = NotAModule("astropy")
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
                log = NotAModule("astropy")
            self._log = log
        return self._log

    _units = None
    @property
    def units(self):
        if self._units is None:
            try:
                from astropy import units
            except ImportError:
                units = None
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
                conv = None
            self._conv = conv
        return self._conv

_astropy = astropy_imports()
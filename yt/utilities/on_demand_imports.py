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

class astropy_imports:
    _pyfits = None
    @property
    def pyfits(self):
        if self._pyfits is None:
            try:
                import astropy.io.fits as pyfits
                self.log
            except ImportError:
                pyfits = None
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
                pywcs = None
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
                log = None
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
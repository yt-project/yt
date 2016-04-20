"""
A wrapper class for h5py file objects.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
from distutils.version import LooseVersion

class HDF5FileHandler(object):
    handle = None

    def __init__(self, filename):
        self.handle = h5py.File(filename, 'r')

    def __del__(self):
        if self.handle is not None:
            # In h5py 2.4 and newer file handles are closed automatically.
            # so only close the handle on old versions.  This works around an
            # issue in h5py when run under python 3.4.
            # See https://github.com/h5py/h5py/issues/534
            if LooseVersion(h5py.__version__) < '2.4.0':
                self.handle.close()

    def __getitem__(self, key):
        return self.handle[key]

    def __contains__(self, item):
        return item in self.handle

    def __len__(self):
        return len(self.handle)

    @property
    def attrs(self):
        return self.handle.attrs

    def keys(self):
        return list(self.handle.keys())

    def items(self):
        return list(self.handle.items())

    def close(self):
        if self.handle is not None:
            self.handle.close()

class FITSFileHandler(HDF5FileHandler):
    def __init__(self, filename):
        from yt.utilities.on_demand_imports import _astropy
        if isinstance(filename, _astropy.pyfits.hdu.image._ImageBaseHDU):
            self.handle = _astropy.pyfits.HDUList(filename)
        elif isinstance(filename, _astropy.pyfits.HDUList):
            self.handle = filename
        else:
            self.handle = _astropy.pyfits.open(
                filename, memmap=True, do_not_scale_image_data=True,
                ignore_blank=True)
        self._fits_files = []

    def __del__(self):
        for f in self._fits_files:
            f.close()
        del self._fits_files
        del self.handle
        self.handle = None
        super(FITSFileHandler, self).__del__()

    def close(self):
        self.handle.close()

class NetCDF4FileHandler(object):
    def __init__(self, filename):
        from netCDF4 import Dataset
        ds = Dataset(filename)
        self.dataset = ds

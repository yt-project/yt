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

import h5py

class HDF5FileHandler(object):
    def __init__(self, filename):
        self.handle = h5py.File(filename, 'r')

    def __del__(self):
        if hasattr(self, 'handle'):
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

    @property
    def keys(self):
        return self.handle.keys

    @property
    def items(self):
        return self.handle.items

class FITSFileHandler(HDF5FileHandler):
    def __init__(self, filename):
        from yt.utilities.on_demand_imports import _astropy
        if isinstance(filename, _astropy.pyfits.PrimaryHDU):
            self.handle = _astropy.pyfits.HDUList(filename)
        else:
            self.handle = _astropy.pyfits.open(
                filename, memmap=True, do_not_scale_image_data=True,
                ignore_blank=True)

    def __del__(self):
        for f in self._fits_files:
            f.close()
            del f
        super(FITSFileHandler, self).__del__()

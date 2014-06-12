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

class FileHandler(object):
    def __init__(self, filename):
        self.handle = h5py.File(filename, 'r')

    def __del__(self):
        self.handle.close()


"""
A set of utilities for analyzing disks

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.mods import *

class StackedDiskImage(object):
    def __init__(self, pf, center, norm_vec, thickness, width,
                 nslices = 100, buff_size=(800,800)):
        self.data = {}
        self.pf = pf
        self.center = na.array(center)
        self.norm_vec = norm_vec/na.sqrt(na.dot(norm_vec,norm_vec))
        self.thickness = thickness
        self.width = width
        self.nslices = nslices
        self.buff_size = buff_size

    def __getitem__(self, key):
        if key not in self.data:
            self[key] = self.get_data(key)
        return self.data[key]

    def __delitem__(self, key):
        del self.data[key]

    def __setitem__(self, key, val):
        self.data[key] = val

    def __iter__(self):
        kk = self.data.keys()
        for key in kk: yield key

    def get_data(self, field):
        # Needs to be rewritten to get multiple fields
        buffer = na.zeros(self.buff_size, dtype='float64')
        wsp = (-self.width/2.0, self.width/2.0,
               -self.width/2.0, self.width/2.0)
        for d in na.mgrid[-self.thickness/2.0:self.thickness/2.0:self.nslices*1j]:
            print d, self.thickness/2.0
            center = self.center + self.norm_vec*d
            cut = self.pf.h.cutting(self.norm_vec, center)
            buff = raven.ObliqueFixedResolutionBuffer(
                        cut, wsp, self.buff_size)
            buffer += buff[field]
        return buffer

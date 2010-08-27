"""
A plotting mechanism based on the idea of a "window" into the data.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from image_writer import \
    write_image
from fixed_resolution import \
    FixedResolutionBuffer

class PlotWindow(object):
    def __init__(self, data_source, bounds, buff_size=(800,800), antialias = True):
        r"""
        PlotWindow(data_source, bounds, buff_size=(800,800), antialias = True)
        
        A ploting mechanism based around the concept of a window into a
        data source. It can have arbitrary fields, each of which will be
        centered on the same viewpoint, but will have individual zlimits. 
        
        The data and plot are updated separately, and each can be
        invalidated as the object is modified.
        
        Data is handled by a FixedResolutionBuffer object.
        """

        self._plot_valid = False
        self._data_valid = False

        self._frb = FixedResolutionBuffer(data_source, 
                                          bounds, buff_size, antialias)
        self._frb._get_data_source_fields()

    def __getitem__(self, item):
        self._frb[item]

    @property
    def fields(self):
        return self._frb.fields

    def save(self):
        """
        REPLACE THIS
        """
        for field in self._frb.data.keys():
            name = "%s.png" % field
            print "writing %s" % name
            write_image(self._frb[field],name)

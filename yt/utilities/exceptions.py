"""
This is a library of yt-defined exceptions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

# We don't need to import 'exceptions'
#import exceptions

class YTException(Exception):
    def __init__(self, pf = None):
        Exception.__init__(self)
        self.pf = pf

# Data access exceptions:

class YTSphereTooSmall(YTException):
    def __init__(self, pf, radius, smallest_cell):
        YTException.__init__(self, pf)
        self.radius = radius
        self.smallest_cell = smallest_cell

    def __str__(self):
        return "%0.5e < %0.5e" % (self.radius, self.smallest_cell)

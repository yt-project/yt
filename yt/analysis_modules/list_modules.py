"""
A mechanism for listing available analysis modules.

Author: Matthew Turk <matthewturk@gmail.com>
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

import os
import sys

def get_available_modules():
    modpath = os.path.abspath(os.path.dirname(__file__))
    available_modules = []
    for d in [os.path.join(modpath, f) for f in os.listdir(modpath)]:
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, "api.py")):
            available_modules.append(os.path.basename(d))
    return available_modules

class AnalysisModuleLoader(object):

    @property
    def available_modules(self):
        return get_available_modules()

    def __getattr__(self, attr):
        try:
            name = "yt.analysis_modules.%s.api" % (attr)
            nm = __import__(name, level=-1)
            setattr(self, attr, sys.modules[name])
        except ImportError:
            raise AttributeError(attr)
        return getattr(self, attr)

    def __dir__(self):
        # This is a badly behaving object.  I was unable to get this line:
        #return super(AnalysisModuleLoader, self).__dir__() + self.available_modules
        # to work, so we simply return only the methods we know about.
        return ["available_modules"] + self.available_modules

amods = AnalysisModuleLoader()

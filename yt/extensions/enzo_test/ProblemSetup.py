"""
This is a definition of a class for defining problem verification.

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

from yt.mods import *
from yt.extensions.EnzoSimulation import *
from ProblemVerification import verification_registry

from mercurial import commands, hg, ui

class Problem(object):
    def __init__(self, name, pfname, repo_url=None):
        self.name = name
        self.pfname = pfname
        self.simulation = EnzoSimulation(pfname)
        if repo_url is None: repo_url = "results_%s" % name
        self.tests = []
        self.init_repo(repo_url)

    def init_repo(self, url):
        pass
        
    def add_test(self, ptype, name, *args, **kwargs):
        pass
        

# repo['tip']['yt/lagos/__init__.py'].data()

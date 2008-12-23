"""
A simple SQLite interface to grabbing and storing parameter files

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

## This is going to be the means of interacting with an SQLite (no, I repeat,
## ORM will be used here!) db in the ~/.yt/ directory where parameter files
## will be stored.  It will be queried for a hash, which then gets instantiated
## and returned to the user.

## Non-functional right now.

from yt.fido import *
import sqlite

class ParameterFileStore(object):
    _shared_state = {} # Do we need this?

    def __new__(cls, *args, **kwargs):
        pass
        
    def _get_db_name(self):
        return os.path.expanduser("~/.yt/pfdb.sql")

    def _initialize_new(self, filename = None):
        self._conn = sqlite.connect(self._get_db_name())

    def get_pf(self, hash):
        pass

    def check_pf(self, fn = None, ctid = None, basedir = None,
                 initial_time = None, hash = None):
        pass

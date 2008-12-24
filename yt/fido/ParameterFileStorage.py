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

## Our table layout:
##  ParameterFiles
##      Filename        fn      text
##      Last known path path    text
##      Sim time        time    real    
##      CurrentTimeID   ctid    real
##      Hash            hash    text

from yt.fido import *
import sqlite3
import yt.lagos.OutputTypes
import os.path

#sqlite3.register_adapter(yt.lagos.OutputTypes.EnzoStaticOutput, _adapt_pf)
#sqlite3.register_converter("pfile", _convert_pf)

class ParameterFileStore(object):

    _shared_state = {}
    _conn = None

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self, in_memory = False):
        if self._conn is None:
            self._conn = sqlite3.connect(self._get_db_name(),
                    detect_types=sqlite3.PARSE_DECLTYPES)
            self._conn.row_factory = self._convert_pf
            self._initialize_new()
        
    def _get_db_name(self):
        return os.path.expanduser("~/.yt/pfdb.sql")

    def _initialize_new(self, filename = None):
        c = self._conn.cursor()
        try:
            c.execute("""create table parameter_files
                            (pf text, path text, time real,
                             ctid real, hash text primary key unique)""")
            self._conn.commit()
        except sqlite3.OperationalError:
            pass
        c.close()

    def wipe_hash(self, hash):
        c = self._conn.cursor()
        c.execute("delete from parameter_files where hash=?", (hash,))
        self._conn.commit()
        c.close()

    def get_pf_hash(self, hash):
        c = self._conn.cursor()
        c.execute("""select * from parameter_files where hash=?""",
                  (hash,))
        return c.fetchall()[0] # Unique

    def get_pf_ctid(self, ctid):
        c = self._conn.cursor()
        c.execute("""select * from parameter_files where ctid=?""",
                  (ctid,))
        return c.fetchall()

    def get_count_hash(self, hash):
        c = self._conn.cursor()
        c.execute("""select pf, path from parameter_files where hash=?""",
                  (hash,))
        res = c.fetchall()
        return len(res), res

    def _adapt_pf(self, pf):
        return (pf.basename, pf.fullpath,
                pf["InitialTime"], pf["CurrentTimeIdentifier"],
                pf._hash())

    def _convert_pf(self, cursor, row):
        if len(row) != 5: return row
        bn, fp, t1, ctid, hash = row
        fn = os.path.join(fp, bn)
        if os.path.exists(fn):
            pf = yt.lagos.EnzoStaticOutput(
                os.path.join(fp, bn))
        else:
            raise IOError
        return pf

    def check_pf(self, pf):
        rc, res = self.get_count_hash(pf._hash())
        if rc == 0:
            self.insert_pf(pf)
            return
        elif rc > 1:
            self.wipe_hash(pf._hash())
            self.insert_pf(pf)
            return
        bn, fp = res[0]
        if bn != pf.basename or fp != pf.fullpath:
            self.wipe_hash(pf._hash())
            self.insert_pf(pf)

    def insert_pf(self, pf):
        c = self._conn.cursor()
        c.execute("""insert into parameter_files values
                     (?,?,?,?,?)""", self._adapt_pf(pf))
        self._conn.commit()
        c.close()

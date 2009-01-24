"""
A simple SQLite interface to grabbing and storing parameter files

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

from yt.config import ytcfg
from yt.fido import *
from yt.funcs import *
import shelve
import os.path

class NoParameterShelf(Exception):
    pass

class ParameterFileStore(object):

    _shared_state = {}
    _shelf = None

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self, in_memory = False):
        self.__init_shelf()

    def _get_db_name(self):
        if not os.access(os.path.expanduser("~/"), os.W_OK):
            return os.path.abspath("parameter_files.db")
        return os.path.expanduser("~/.yt/parameter_files.db")

    def wipe_hash(self, hash):
        if hash in self.keys():
            del self[hash]

    def get_pf_hash(self, hash):
        return self._convert_pf(self[hash])

    def get_pf_ctid(self, ctid):
        for h in self.keys():
            if self[h]['ctid'] == ctid:
                return self._convert_pf(self[h])

    def _adapt_pf(self, pf):
        return dict(bn=pf.basename,
                    fp=pf.fullpath,
                    tt=pf["InitialTime"],
                    ctid=pf["CurrentTimeIdentifier"])

    def _convert_pf(self, pf_dict):
        bn = pf_dict['bn']
        fp = pf_dict['fp']
        fn = os.path.join(fp, bn)
        mylog.info("Checking %s", fn)
        if os.path.exists(fn):
            import yt.lagos.OutputTypes as ot
            pf = ot.EnzoStaticOutput(
                os.path.join(fp, bn))
        else:
            raise IOError
        return pf

    def check_pf(self, pf):
        if pf._hash() not in self.keys():
            self.insert_pf(pf)
            return
        pf_dict = self[pf._hash()]
        if pf_dict['bn'] != pf.basename \
          or pf_dict['fp'] != pf.fullpath:
            self.wipe_hash(pf._hash())
            self.insert_pf(pf)

    def __read_only(self):
        if self._shelf is not None: return self._shelf
        return shelve.open(self._get_db_name(), flag='r', protocol=-1)

    def __read_write(self):
        if self._shelf is not None: return self._shelf
        return shelve.open(self._get_db_name(), flag='c', protocol=-1)

    def insert_pf(self, pf):
        self[pf._hash()] = self._adapt_pf(pf)

    def __getitem__(self, key):
        my_shelf = self.__read_only()
        return my_shelf[key]

    def __store_item(self, key, val):
        my_shelf = self.__read_write()
        my_shelf[key] = val

    def __delete_item(self, key):
        my_shelf = self.__read_write()
        del my_shelf[key]

    def __init_shelf(self):
        dbn = self._get_db_name()
        dbdir = os.path.dirname(dbn)
        if not ytcfg.getboolean("yt", "StoreParameterFiles"):
            # This ensures that even if we're not storing them in the file
            # system, we're at least keeping track of what we load
            self._shelf = defaultdict(lambda: dict(bn='',fp='',tt='',ctid=''))
            return
        try:
            if not os.path.isdir(dbdir): os.mkdir(dbdir)
        except OSError:
            raise NoParameterShelf()
        only_on_root(shelve.open, dbn, 'c', protocol=-1)

    def __setitem__(self, key, val):
        only_on_root(self.__store_item, key, val)

    def __delitem__(self, key):
        only_on_root(self.__delete_item, key)

    def keys(self):
        my_shelf = self.__read_only()
        return my_shelf.keys()

class ObjectStorage(object):
    pass

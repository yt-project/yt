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
from yt.lagos.ParallelTools import parallel_simple_proxy
import csv
import os.path

output_type_registry = {}
_field_names = ('hash','bn','fp','tt','ctid','class_name')

class NoParameterShelf(Exception):
    pass

class UnknownStaticOutputType(Exception):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return "%s" % self.name

    def __repr__(self):
        return "%s" % self.name

class ParameterFileStore(object):

    _shared_state = {}
    _distributed = True
    _processing = False
    _owner = 0

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self, in_memory = False):
        if ytcfg.getboolean("yt", "StoreParameterFiles"):
            self._read_only = False
            self.init_db()
            self._records = self.read_db()
        else:
            self._read_only = True
            self._records = {}

    @parallel_simple_proxy
    def init_db(self):
        dbn = self._get_db_name()
        dbdir = os.path.dirname(dbn)
        try:
            if not os.path.isdir(dbdir): os.mkdir(dbdir)
        except OSError:
            raise NoParameterShelf()
        open(dbn, 'ab') # make sure it exists, allow to close
        # Now we read in all our records and return them
        # these will be broadcast

    def _get_db_name(self):
        base_file_name = ytcfg.get("yt","ParameterFileStore")
        if not os.access(os.path.expanduser("~/"), os.W_OK):
            return os.path.abspath(base_file_name)
        return os.path.expanduser("~/.yt/%s" % base_file_name)

    def get_pf_hash(self, hash):
        return self._convert_pf(self._records[hash])

    def get_pf_ctid(self, ctid):
        for h in self._records:
            if self._records[h]['ctid'] == ctid:
                return self._convert_pf(self._records[h])

    def _adapt_pf(self, pf):
        return dict(bn=pf.basename,
                    fp=pf.fullpath,
                    tt=pf["InitialTime"],
                    ctid=pf["CurrentTimeIdentifier"],
                    class_name=pf.__class__.__name__)

    def _convert_pf(self, pf_dict):
        bn = pf_dict['bn']
        fp = pf_dict['fp']
        fn = os.path.join(fp, bn)
        class_name = pf_dict['class_name']
        if class_name not in output_type_registry:
            raise UnknownStaticOutputType(class_name)
        mylog.info("Checking %s", fn)
        if os.path.exists(fn):
            pf = output_type_registry[class_name](os.path.join(fp, bn))
        else:
            raise IOError
        return pf

    def check_pf(self, pf):
        if pf._hash() not in self._records:
            self.insert_pf(pf)
            return
        pf_dict = self._records[pf._hash()]
        if pf_dict['bn'] != pf.basename \
          or pf_dict['fp'] != pf.fullpath:
            self.wipe_hash(pf._hash())
            self.insert_pf(pf)

    def insert_pf(self, pf):
        self._records[pf._hash()] = self._adapt_pf(pf)
        self.flush_db()

    def wipe_hash(self, hash):
        if hash not in self._records: return
        del self._records[hash]
        self.flush_db()

    def flush_db(self):
        if self._read_only: return
        self._write_out()
        self.read_db()

    @parallel_simple_proxy
    def _write_out(self):
        if self._read_only: return
        f = open(self._get_db_name(), 'ab')
        f.seek(0,2)
        w = csv.DictWriter(f, _field_names)
        for h,v in sorted(self._records.items()):
            v['hash'] = h
            w.writerow(v)
        f.close()

    @parallel_simple_proxy
    def read_db(self):
        f=open(self._get_db_name(), 'rb')
        vals = csv.DictReader(f, _field_names)
        db = {}
        for v in vals:
            db[v.pop('hash')] = v
            
        return db

class ObjectStorage(object):
    pass

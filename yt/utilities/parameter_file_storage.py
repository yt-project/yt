"""
A simple CSV database for grabbing and storing datasets



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import csv
import os.path
from itertools import islice

from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_simple_proxy

output_type_registry = {}
simulation_time_series_registry = {}
_field_names = ('hash', 'bn', 'fp', 'tt', 'ctid', 'class_name', 'last_seen')

class NoParameterShelf(Exception):
    pass

class UnknownDatasetType(Exception):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return "%s" % self.name

    def __repr__(self):
        return "%s" % self.name

class ParameterFileStore(object):
    """
    This class is designed to be a semi-persistent storage for parameter
    files.  By identifying each dataset with a unique hash, objects
    can be stored independently of datasets -- when an object is
    loaded, the dataset is as well, based on the hash.  For
    storage concerns, only a few hundred will be retained in cache.

    """

    _shared_state = {}
    _distributed = True
    _processing = False
    _owner = 0
    _register = True

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self, in_memory=False):
        """
        Create the dataset database if yt is configured to store them.
        Otherwise, use read-only settings.

        """
        if self._register is False: return
        if ytcfg.getboolean("yt", "StoreParameterFiles"):
            self._read_only = False
            self.init_db()
            self._records = self.read_db()
        else:
            self._read_only = True
            self._records = {}
        self._register = False

    @parallel_simple_proxy
    def init_db(self):
        """
        This function ensures that the storage database exists and can be used.
        """
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
        base_file_name = ytcfg.get("yt", "ParameterFileStore")
        if not os.access(os.path.expanduser("~/"), os.W_OK):
            return os.path.abspath(base_file_name)
        return os.path.expanduser("~/.yt/%s" % base_file_name)

    def get_ds_hash(self, hash):
        """ This returns a dataset based on a hash. """
        return self._convert_ds(self._records[hash])

    def get_ds_ctid(self, ctid):
        """ This returns a dataset based on a CurrentTimeIdentifier. """
        for h in self._records:
            if self._records[h]['ctid'] == ctid:
                return self._convert_ds(self._records[h])

    def _adapt_ds(self, ds):
        """ This turns a dataset into a CSV entry. """
        return dict(bn=ds.basename,
                    fp=ds.fullpath,
                    tt=ds.current_time,
                    ctid=ds.unique_identifier,
                    class_name=ds.__class__.__name__,
                    last_seen=ds._instantiated)

    def _convert_ds(self, ds_dict):
        """ This turns a CSV entry into a dataset. """
        bn = ds_dict['bn']
        fp = ds_dict['fp']
        fn = os.path.join(fp, bn)
        class_name = ds_dict['class_name']
        if class_name not in output_type_registry:
            raise UnknownDatasetType(class_name)
        mylog.info("Checking %s", fn)
        if os.path.exists(fn):
            ds = output_type_registry[class_name](os.path.join(fp, bn))
        else:
            raise IOError
        # This next one is to ensure that we manually update the last_seen
        # record *now*, for during write_out.
        self._records[ds._hash()]['last_seen'] = ds._instantiated
        return ds

    def check_ds(self, ds):
        """
        This will ensure that the dataset (*ds*) handed to it is
        recorded in the storage unit.  In doing so, it will update path
        and "last_seen" information.
        """
        hash = ds._hash()
        if hash not in self._records:
            self.insert_ds(ds)
            return
        ds_dict = self._records[hash]
        self._records[hash]['last_seen'] = ds._instantiated
        if ds_dict['bn'] != ds.basename \
          or ds_dict['fp'] != ds.fullpath:
            self.wipe_hash(hash)
            self.insert_ds(ds)

    def insert_ds(self, ds):
        """ This will insert a new *ds* and flush the database to disk. """
        self._records[ds._hash()] = self._adapt_ds(ds)
        self.flush_db()

    def wipe_hash(self, hash):
        """
        This removes a *hash* corresponding to a dataset from the
        storage.
        """
        if hash not in self._records: return
        del self._records[hash]
        self.flush_db()

    def flush_db(self):
        """ This flushes the storage to disk. """
        if self._read_only: return
        self._write_out()
        self.read_db()

    def get_recent(self, n=10):
        recs = sorted(self._records.values(), key=lambda a: -a['last_seen'])[:n]
        return recs

    @parallel_simple_proxy
    def _write_out(self):
        if self._read_only: return
        fn = self._get_db_name()
        f = open("%s.tmp" % fn, 'wb')
        w = csv.DictWriter(f, _field_names)
        maxn = ytcfg.getint("yt","maximumstoreddatasets") # number written
        for h,v in islice(sorted(self._records.items(),
                          key=lambda a: -a[1]['last_seen']), 0, maxn):
            v['hash'] = h
            w.writerow(v)
        f.close()
        os.rename("%s.tmp" % fn, fn)

    @parallel_simple_proxy
    def read_db(self):
        """ This will read the storage device from disk. """
        f = open(self._get_db_name(), 'rb')
        vals = csv.DictReader(f, _field_names)
        db = {}
        for v in vals:
            db[v.pop('hash')] = v
            if v['last_seen'] is None:
                v['last_seen'] = 0.0
            else: v['last_seen'] = float(v['last_seen'])
        return db

class ObjectStorage(object):
    pass

class EnzoRunDatabase(object):
    conn = None

    def __init__(self, path = None):
        if path is None:
            path = ytcfg.get("yt", "enzo_db")
            if len(path) == 0:
                raise RuntimeError
        import sqlite3
        self.conn = sqlite3.connect(path)

    def find_uuid(self, u):
        cursor = self.conn.execute(
            "select ds_path from enzo_outputs where dset_uuid = '%s'" % (
                u))
        # It's a 'unique key'
        result = cursor.fetchone()
        if result is None: return None
        return result[0]

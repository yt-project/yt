"""
A simple CSV database for grabbing and storing parameter files

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

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

import csv
import os.path
from itertools import islice

from yt.config import ytcfg
from yt.funcs import *
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_simple_proxy

import yt.utilities.peewee as peewee

output_type_registry = {}
_field_names = ('hash', 'bn', 'fp', 'tt', 'ctid', 'class_name', 'last_seen')

class NoParameterShelf(Exception):
    pass

class UnknownStaticOutputType(Exception):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return "%s" % self.name

    def __repr__(self):
        return "%s" % self.name

_field_spec = dict(
    dset_uuid = peewee.TextField(),
    output_type = peewee.TextField(),
    pf_path = peewee.TextField(),
    creation_time = peewee.IntegerField(),
    last_seen_time = peewee.IntegerField(),
    simulation_uuid = peewee.TextField(),
    redshift = peewee.FloatField(),
    time = peewee.FloatField(),
    topgrid0 = peewee.IntegerField(),
    topgrid1 = peewee.IntegerField(),
    topgrid2 = peewee.IntegerField(),
)

class ParameterFileStore(object):
    """
    This class is designed to be a semi-persistent storage for parameter
    files.  By identifying each parameter file with a unique hash, objects
    can be stored independently of parameter files -- when an object is
    loaded, the parameter file is as well, based on the hash.  For
    storage concerns, only a few hundred will be retained in cache.

    """

    _shared_state = {}
    _distributed = True
    _processing = False
    _owner = 0
    conn = None

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self, in_memory=False):
        """
        Create the parameter file database if yt is configured to store them.
        Otherwise, use read-only settings.

        """
        if ytcfg.getboolean("yt", "StoreParameterFiles"):
            self._read_only = False
            self.init_db()
        else:
            self._read_only = True
            self._records = {}

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
        self.conn = peewee.SqliteDatabase(dbn)
        class SimulationOutputsMeta:
            database = self.conn
            db_table = "simulation_outputs"
        _field_spec["Meta"] = SimulationOutputsMeta
        self.output_model = type(
            "SimulationOutputs",
            (peewee.Model,),
            _field_spec,
        )
        self.output_model._meta.pk_name = "dset_uuid"
        try:
            self.conn.connect()
        except:
            self.conn = None
        try:
            self.output_model.create_table()
        except:
            pass
        self.conn = None

    def _get_db_name(self):
        base_file_name = ytcfg.get("yt", "ParameterFileStore")
        if not os.access(os.path.expanduser("~/"), os.W_OK):
            return os.path.abspath(base_file_name)
        return os.path.expanduser("~/.yt/%s" % base_file_name)

    def get_pf_hash(self, hash):
        if self.conn is None: return
        """ This returns a parameter file based on a hash. """
        output = self.output_model.get(dset_uuid = hash)
        return self._convert_pf(output)

    def _convert_pf(self, inst):
        """ This turns a model into a parameter file. """
        if self.conn is None: return
        fn = inst.pf_path
        if inst.output_type not in output_type_registry:
            raise UnknownStaticOutputType(inst.output_type)
        mylog.info("Checking %s", fn)
        if os.path.exists(fn):
            pf = output_type_registry[inst.output_type](fn)
        else:
            raise IOError
        # This next one is to ensure that we manually update the last_seen
        # record *now*, for during write_out.
        self.output_model.update(last_seen_time = pf._instantiated).where(
            dset_uuid = inst.dset_uuid).execute()
        return pf

    def check_pf(self, pf):
        """
        This will ensure that the parameter file (*pf*) handed to it is
        recorded in the storage unit.  In doing so, it will update path
        and "last_seen" information.
        """
        if self.conn is None: return
        q = self.output_model.select().where(dset_uuid = pf._hash())
        q.execute()
        if q.count() == 0:
            self.insert_pf(pf)
            return
        # Otherwise we update
        self.output_model.update(
            last_seen_time = pf._instantiated,
            pf_path = os.path.join(pf.basename, pf.fullpath)
        ).where(
            dset_uuid = pf._hash()).execute(
        )

    def insert_pf(self, pf):
        """ This will insert a new *pf* and flush the database to disk. """
        if self.conn is None: return
        q = self.output_model.insert(
                    dset_uuid = pf._hash(),
                    output_type = pf.__class__.__name__,
                    pf_path = os.path.join(
                        pf.fullpath, pf.basename),
                    creation_time = pf.parameters.get(
                        "CurrentTimeIdentifier", 0), # Get os.stat
                    last_seen_time = pf._instantiated,
                    simulation_uuid = pf.parameters.get(
                        "SimulationUUID", ""), # NULL
                    redshift = pf.current_redshift,
                    time = pf.current_time,
                    topgrid0 = pf.domain_dimensions[0],
                    topgrid1 = pf.domain_dimensions[1],
                    topgrid2 = pf.domain_dimensions[2])
        q.execute()

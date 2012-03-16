"""
Operations to get Rockstar loaded up

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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
from os import environ
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, ProcessorPool, Communicator

import rockstar_interface
import socket
import time

class DomainDecomposer(ParallelAnalysisInterface):
    def __init__(self, pf, comm):
        ParallelAnalysisInterface.__init__(self, comm=comm)
        self.pf = pf
        self.hierarchy = pf.h
        self.center = (pf.domain_left_edge + pf.domain_right_edge)/2.0

    def decompose(self):
        dd = self.pf.h.all_data()
        check, LE, RE, data_source = self.partition_hierarchy_3d(dd)
        return data_source

class RockstarHaloFinder(ParallelAnalysisInterface):
    def __init__(self, pf, num_readers = 0, num_writers = 0):
        ParallelAnalysisInterface.__init__(self)
        # No subvolume support
        self.pf = pf
        self.hierarchy = pf.h
        self.num_readers = num_readers
        self.num_writers = num_writers
        if self.num_readers + self.num_writers + 1 != self.comm.size:
            raise RuntimeError
        self.center = (pf.domain_right_edge + pf.domain_left_edge)/2.0
        data_source = None
        if self.comm.size > 1:
            self.pool = ProcessorPool()
            self.pool.add_workgroup(1, name = "server")
            self.pool.add_workgroup(num_readers, name = "readers")
            self.pool.add_workgroup(num_writers, name = "writers")
            for wg in self.pool.workgroups:
                if self.comm.rank in wg.ranks: self.workgroup = wg
        data_source = self.pf.h.all_data()
        self.handler = rockstar_interface.RockstarInterface(
                self.pf, data_source)

    def _get_hosts(self):
        if self.comm.size == 1 or self.workgroup.name == "server":
            server_address = socket.gethostname()
            sock = socket.socket()
            sock.bind(('', 0))
            port = sock.getsockname()[-1]
            del sock
        else:
            server_address, port = None, None
        self.server_address, self.port = self.comm.mpi_bcast_pickled(
            (server_address, port))
        self.port = str(self.port)

    def run(self, block_ratio = 1):
        if block_ratio != 1:
            raise NotImplementedError
        self._get_hosts()
        self.handler.setup_rockstar(self.server_address, self.port,
                    parallel = self.comm.size > 1,
                    num_readers = self.num_readers,
                    num_writers = self.num_writers,
                    writing_port = -1,
                    block_ratio = block_ratio)
        if self.comm.size == 1:
            self.handler.call_rockstar()
        else:
            self.comm.barrier()
            if self.workgroup.name == "server":
                self.handler.start_server()
            elif self.workgroup.name == "readers":
                time.sleep(0.5 + self.workgroup.comm.rank/10.0)
                self.handler.start_client()
            elif self.workgroup.name == "writers":
                time.sleep(1.0 + self.workgroup.comm.rank/10.0)
                self.handler.start_client()
        self.comm.barrier()

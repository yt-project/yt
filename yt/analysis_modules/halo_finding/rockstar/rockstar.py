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
from os import mkdir
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, ProcessorPool, Communicator

from yt.analysis_modules.halo_finding.halo_objects import * #Halos & HaloLists
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
    def __init__(self, pf, num_readers = 1, num_writers = None, 
            outbase=None,particle_mass=-1.0,overwrite=False,
            left_edge = None, right_edge = None):
        ParallelAnalysisInterface.__init__(self)
        # No subvolume support
        self.pf = pf
        self.hierarchy = pf.h
        if num_writers is None:
            num_writers = self.comm.size - num_readers -1
        self.num_readers = num_readers
        self.num_writers = num_writers
        self.particle_mass = particle_mass 
        self.overwrite = overwrite
        if left_edge is None:
            left_edge = pf.domain_left_edge
        if right_edge is None:
            right_edge = pf.domain_right_edge
        self.le = left_edge
        self.re = right_edge
        if self.num_readers + self.num_writers + 1 != self.comm.size:
            print '%i reader + %i writers != %i mpi'%\
                    (self.num_readers, self.num_writers, self.comm.size)
            raise RuntimeError
        self.center = (pf.domain_right_edge + pf.domain_left_edge)/2.0
        data_source = self.pf.h.all_data()
        self.handler = rockstar_interface.RockstarInterface(
                self.pf, data_source)
        if outbase is None:
            outbase = str(self.pf)+'_rockstar'
        self.outbase = outbase        

    def _get_hosts(self):
        if self.comm.size == 1 or self.workgroup.name == "server":
            server_address = socket.gethostname()
            sock = socket.socket()
            sock.bind(('', 0))
            port = sock.getsockname()[-1]
            del sock
        else:
            server_address, port = None, None
        self.server_address, self.port = self.comm.mpi_bcast(
            (server_address, port))
        self.port = str(self.port)

    def run(self, block_ratio = 1,**kwargs):
        """
        
        """
        if self.comm.size > 1:
            self.pool = ProcessorPool()
            mylog.debug("Num Writers = %s Num Readers = %s",
                        self.num_writers, self.num_readers)
            self.pool.add_workgroup(1, name = "server")
            self.pool.add_workgroup(self.num_readers, name = "readers")
            self.pool.add_workgroup(self.num_writers, name = "writers")
            for wg in self.pool.workgroups:
                if self.comm.rank in wg.ranks: self.workgroup = wg
        if block_ratio != 1:
            raise NotImplementedError
        self._get_hosts()
        #because rockstar *always* write to exactly the same
        #out_0.list filename we make a directory for it
        #to sit inside so it doesn't get accidentally
        #overwritten 
        if self.workgroup.name == "server":
            if not os.path.exists(self.outbase):
                os.mkdir(self.outbase)
        self.handler.setup_rockstar(self.server_address, self.port,
                    parallel = self.comm.size > 1,
                    num_readers = self.num_readers,
                    num_writers = self.num_writers,
                    writing_port = -1,
                    block_ratio = block_ratio,
                    outbase = self.outbase,
                    particle_mass = float(self.particle_mass),
                    **kwargs)
        if self.comm.size == 1:
            self.handler.call_rockstar()
        else:
            self.comm.barrier()
            if self.workgroup.name == "server":
                self.handler.start_server()
            elif self.workgroup.name == "readers":
                time.sleep(0.1 + self.workgroup.comm.rank/10.0)
                self.handler.start_client()
            elif self.workgroup.name == "writers":
                time.sleep(0.2 + self.workgroup.comm.rank/10.0)
                self.handler.start_client()
            self.pool.free_all()
        self.comm.barrier()
        #quickly rename the out_0.list 
    
    def halo_list(self,file_name='out_0.list'):
        """
        Reads in the out_0.list file and generates RockstarHaloList
        and RockstarHalo objects.
        """
        return RockstarHaloList(self.pf,self.outbase+'/%s'%file_name)

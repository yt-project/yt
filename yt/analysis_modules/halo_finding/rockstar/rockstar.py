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
    def __init__(self, ts, num_readers = 1, num_writers = None, 
            outbase=None,particle_mass=-1.0,dm_type=1):
        r"""Spawns the Rockstar Halo finder, distributes dark matter
        particles and finds halos.

        The halo finder requires dark matter particles of a fixed size.
        Rockstar has three main processes: reader, writer, and the 
        server which coordinates reader/writer processes.

        Parameters
        ----------
        ts   : TimeSeriesData, StaticOutput
            This is the data source containing the DM particles. Because 
            halo IDs may change from one snapshot to the next, the only
            way to keep a consistent halo ID across time is to feed 
            Rockstar a set of snapshots, ie, via TimeSeriesData.
        num_readers: int
            The number of reader can be increased from the default
            of 1 in the event that a single snapshot is split among
            many files. This can help in cases where performance is
            IO-limited. Default is 1.
        num_writers: int
            The number of writers determines the number of processing threads
            as well as the number of threads writing output data.
            The default is set comm.size-num_readers-1.
        outbase: str
            This is where the out*list files that Rockstar makes should be
            placed. Default is str(pf)+'_rockstar'.
        particle_mass: float
            This sets the DM particle mass used in Rockstar.
        dm_type: 1
            In order to exclude stars and other particle types, define
            the dm_type. Default is 1, as Enzo has the DM particle type=1.

        Returns
        -------
        None

        Examples
        --------
        To use the script below you must run it using MPI:
        mpirun -np 3 python test_rockstar.py --parallel

        test_rockstar.py:

        from mpi4py import MPI
        from yt.analysis_modules.halo_finding.rockstar.api import RockstarHaloFinder
        from yt.mods import *
        import sys

        files = glob.glob('/u/cmoody3/data/a*')
        files.sort()
        ts = TimeSeriesData.from_filenames(files)
        pm = 7.81769027e+11
        rh = RockstarHaloFinder(ts, particle_mass=pm)
        rh.run()
        """
        ParallelAnalysisInterface.__init__(self)
        # No subvolume support
        #we assume that all of the snapshots in the time series
        #use the same domain info as the first snapshots
        if not isinstance(ts,TimeSeriesData):
            ts = TimeSeriesData([ts])
        self.ts = ts
        self.dm_type = dm_type
        if self.comm.size > 1: 
            self.comm.barrier()            
        tpf = ts.__iter__().next()
        def _particle_count(field,data):
            return (data["particle_type"]==0).sum()
        add_field("particle_count",function=_particle_count,particle_type=True)
        #d = tpf.h.all_data()
        #total_particles = dd.quantities['TotalQuantity']("particle_count")
        #mylog.info("Found %i halo particles",total_particles)
        self.total_particles = -1
        self.hierarchy = tpf.h
        self.particle_mass = particle_mass 
        self.center = (tpf.domain_right_edge + tpf.domain_left_edge)/2.0
        data_source = tpf.h.all_data()
        if outbase is None:
            outbase = str(tpf)+'_rockstar'
        self.outbase = outbase        
        if num_writers is None:
            num_writers = self.comm.size - num_readers -1
        self.num_readers = num_readers
        self.num_writers = num_writers
        if self.num_readers + self.num_writers + 1 != self.comm.size:
            #we need readers+writers+1 server = comm size        
            raise RuntimeError
        self.center = (tpf.domain_right_edge + tpf.domain_left_edge)/2.0
        data_source = tpf.h.all_data()
        self.comm.barrier()
        self.handler = rockstar_interface.RockstarInterface(
                self.ts, data_source)

    def __del__(self):
        self.pool.free_all()

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
        self.handler.setup_rockstar(self.server_address, self.port,
                    len(self.ts), #self.total_particles, 
                    self.dm_type,
                    parallel = self.comm.size > 1,
                    num_readers = self.num_readers,
                    num_writers = self.num_writers,
                    writing_port = -1,
                    block_ratio = block_ratio,
                    outbase = self.outbase,
                    particle_mass = float(self.particle_mass),
                    **kwargs)
        #because rockstar *always* write to exactly the same
        #out_0.list filename we make a directory for it
        #to sit inside so it doesn't get accidentally
        #overwritten 
        if self.workgroup.name == "server":
            if not os.path.exists(self.outbase):
                os.mkdir(self.outbase)
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
        self.pool.free_all()
    
    def halo_list(self,file_name='out_0.list'):
        """
        Reads in the out_0.list file and generates RockstarHaloList
        and RockstarHalo objects.
        """
        tpf = self.ts[0]
        return RockstarHaloList(tpf,file_name)

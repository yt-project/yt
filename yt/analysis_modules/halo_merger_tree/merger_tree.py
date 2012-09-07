"""
MergerTree class and member functions.

Author: Stephen Skory <s@skory.us>
Affiliation: CASS/UC San Diego, CA
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Stephen Skory.  All Rights Reserved.

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

import numpy as na
import os, glob, time, gc, md5, sys
import h5py
import types

from yt.funcs import *

from yt.analysis_modules.halo_finding.halo_objects import \
    FOFHaloFinder, HaloFinder, parallelHF
from yt.analysis_modules.halo_profiler.multi_halo_profiler import \
    HaloProfiler
from yt.convenience import load
from yt.utilities.logger import ytLogger as mylog
import yt.utilities.pydot as pydot
from yt.utilities.spatial import cKDTree
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelDummy, \
    ParallelAnalysisInterface, \
    parallel_blocking_call

try:
    import sqlite3 as sql
except ImportError:
    mylog.error("sqlite3 not imported!")

column_types = {
"GlobalHaloID":"INTEGER",
"SnapCurrentTimeIdentifier":"INTEGER",
"SnapZ":"FLOAT",
"SnapHaloID":"INTEGER",
"HaloMass":"FLOAT",
"NumPart":"INTEGER",
"CenMassX":"FLOAT",
"CenMassY":"FLOAT",
"CenMassZ":"FLOAT",
"BulkVelX":"FLOAT",
"BulkVelY":"FLOAT",
"BulkVelZ":"FLOAT",
"MaxRad":"FLOAT",
"ChildHaloID0":"INTEGER",
"ChildHaloFrac0":"FLOAT",
"ChildHaloID1":"INTEGER",
"ChildHaloFrac1":"FLOAT",
"ChildHaloID2":"INTEGER",
"ChildHaloFrac2":"FLOAT",
"ChildHaloID3":"INTEGER",
"ChildHaloFrac3":"FLOAT",
"ChildHaloID4":"INTEGER", 
"ChildHaloFrac4":"FLOAT"}

# In order.
columns = ["GlobalHaloID", "SnapCurrentTimeIdentifier", "SnapZ", 
"SnapHaloID", "HaloMass", "NumPart", "CenMassX", "CenMassY",
"CenMassZ", "BulkVelX", "BulkVelY", "BulkVelZ", "MaxRad",
"ChildHaloID0", "ChildHaloFrac0",
"ChildHaloID1", "ChildHaloFrac1",
"ChildHaloID2", "ChildHaloFrac2",
"ChildHaloID3", "ChildHaloFrac3",
"ChildHaloID4", "ChildHaloFrac4"]

# Below we make the SQL command that creates the table "Halos" in the
# database. This table is where all the data is stored.
# Each column of data is named and its datatype is specified.
# The GlobalHaloID is given the PRIMARY KEY property, which means that
# the SQLite machinery assigns a consecutive and unique integer value
# to that field automatically as each new entry is entered (that is,
# if GlobalHaloID isn't specified already).
create_db_line = "CREATE TABLE Halos ("
for i, col in enumerate(columns):
    if i == 0:
        create_db_line += "%s %s PRIMARY KEY," % (col, column_types[col])
    else:
        create_db_line += " %s %s," % (col, column_types[col])
# Clean of trailing comma, and closing stuff.
create_db_line = create_db_line[:-1] + ");"

NumNeighbors = 15
NumDB = 5

def minus_one():
    return -1

class DatabaseFunctions(object):
    # Common database functions so it doesn't have to be repeated.
    def _open_database(self):
        # open the database. Check to make sure the database file exists.
        if not os.path.exists(self.database):
            mylog.error("The database file %s cannot be found. Exiting." % \
                self.database)
            return False
        self.conn = sql.connect(self.database)
        self.cursor = self.conn.cursor()
        return True

    def _close_database(self):
        # close the database cleanly.
        self.cursor.close()
        self.conn.close()

class MergerTree(DatabaseFunctions, ParallelAnalysisInterface):
    def __init__(self, restart_files=[], database='halos.db',
            halo_finder_function=HaloFinder, halo_finder_threshold=80.0,
            FOF_link_length=0.2, dm_only=False, refresh=False,
            index=True):
        r"""Build a merger tree of halos over a time-ordered set of snapshots.
        This will run a halo finder to find the halos first if it hasn't already
        been done. The output is a SQLite database file, which may need to
        be stored on a different disk than the data snapshots. See the full
        documentation for details.
        
        Parameters
        ---------
        restart_files : List of strings
            A list containing the paths to the forward time-ordered set of
            data snapshots.
        database : String
            Name of SQLite database file. Default = "halos.db".
        halo_finder_function : HaloFinder name
            The name of the halo finder to use if halo finding is run by 
            the merger tree. Options: HaloFinder, FOFHaloFinder, parallelHF.
            Note that this is not a string, so no quotes. Default = HaloFinder.
        halo_finder_threshold : Float
            If using HaloFinder or parallelHF, the value of the density threshold
            used when halo finding. Default = 80.0.
        FOF_link_length : Float
            If using FOFHaloFinder, the linking length between particles.
            Default = 0.2.
        dm_only : Boolean
            When halo finding, whether to restrict to only dark matter particles.
            Default = False.
        refresh : Boolean
            True forces the halo finder to run even if the halo data has been
            detected on disk. Default = False.
        index : Boolean
            SQLite databases can have added to them an index which greatly
            speeds up future queries of the database,
            at the cost of doubling the disk space used by the file.
            Default = True.

        Examples:
        >>> rf = ['/scratch/user/sim1/DD0000/data0000',
        ... '/scratch/user/sim1/DD0001/data0001',
        ... '/scratch/user/sim1/DD0002/data0002']
        >>> MergerTree(rf, database = '/home/user/sim1-halos.db',
        ... halo_finder_function=parallelHF)
        """
        ParallelAnalysisInterface.__init__(self)
        self.restart_files = restart_files # list of enzo restart files
        self.with_halos = na.ones(len(restart_files), dtype='bool')
        self.database = database # the sqlite database of haloes.
        self.halo_finder_function = halo_finder_function # which halo finder to use
        self.halo_finder_threshold = halo_finder_threshold # overdensity threshold
        self.FOF_link_length= FOF_link_length # For FOF
        self.dm_only = dm_only
        self.refresh = refresh
        self.index = index
        self.zs = {}
        # MPI stuff
        if self.comm.rank is None:
            self.comm.rank = 0
        if self.comm.size is None:
            self.comm.size = 1
        # Get to work.
        if self.refresh and self.comm.rank == 0:
            try:
                os.unlink(self.database)
            except:
                pass
        if self.comm.rank == 0:
            self._open_create_database()
            self._create_halo_table()
        self._run_halo_finder_add_to_db()
        # Find the h5 file names for all the halos.
        for snap in self.restart_files:
            self._build_h5_refs(snap)
        # Find out how much work is already stored in the database.
        if self.comm.rank == 0:
            z_progress = self._find_progress()
        else:
            z_progress = None
        z_progress = self.comm.mpi_bcast(z_progress)
        # Loop over the pairs of snapshots to locate likely neighbors, and
        # then use those likely neighbors to compute fractional contributions.
        last = None
        self.write_values = []
        self.write_values_dict = defaultdict(dict)
        for snap, pair in enumerate(zip(self.restart_files[:-1], self.restart_files[1:])):
            if not self.with_halos[snap] or not self.with_halos[snap+1]:
                continue
            if self.zs[pair[0]] > z_progress:
                continue
            self._find_likely_children(pair[0], pair[1])
            # last is the data for the parent dataset, which can be supplied
            # as the child from the previous round for all but the first loop.
            last = self._compute_child_fraction(pair[0], pair[1], last)
            if self.comm.rank == 0:
                mylog.info("Updating database with parent-child relationships.")
                self._copy_and_update_db()
                # This has to happen because we delete the old database above.
                self._open_create_database()
        del last
        if self.comm.rank == 0:
            if self.index:
                self._write_index()
            self._close_database()
        self.comm.barrier()
        mylog.info("Done!")
        
    def _read_halo_lists(self):
        self.halo_lists = []
        for i,file in enumerate(self.halo_files):
            hp = HaloProfiler(self.restart_files[i], halo_list_file=file)
            self.halo_lists.append(hp.all_halos)

    def _run_halo_finder_add_to_db(self):
        for cycle, file in enumerate(self.restart_files):
            gc.collect()
            pf = load(file)
            self.zs[file] = pf.current_redshift
            self.period = pf.domain_right_edge - pf.domain_left_edge
            # If the halos are already found, skip this data step, unless
            # refresh is True.
            dir = os.path.dirname(file)
            if os.path.exists(os.path.join(dir, 'MergerHalos.out')) and \
                    os.path.exists(os.path.join(dir, 'MergerHalos.txt')) and \
                    glob.glob(os.path.join(dir, 'MergerHalos*h5')) is not [] and \
                    not self.refresh:
                pass
            else:
                # Run the halo finder.
                if self.halo_finder_function == FOFHaloFinder:
                    halos = self.halo_finder_function(pf,
                        link=self.FOF_link_length, dm_only=self.dm_only)
                else:
                    halos = self.halo_finder_function(pf,
                        threshold=self.halo_finder_threshold, dm_only=self.dm_only)
                halos.write_out(os.path.join(dir, 'MergerHalos.out'))
                halos.write_particle_lists(os.path.join(dir, 'MergerHalos'))
                halos.write_particle_lists_txt(os.path.join(dir, 'MergerHalos'))
                if len(halos) == 0:
                    mylog.info("Dataset %s has no halos." % file)
                    self.with_halos[cycle] = False
                    continue
                del halos
            # Now add halo data to the db if it isn't already there by
            # checking the first halo.
            continue_check = False
            if self.comm.rank == 0:
                currt = pf.unique_identifier
                line = "SELECT GlobalHaloID from Halos where SnapHaloID=0\
                and SnapCurrentTimeIdentifier=%d;" % currt
                self.cursor.execute(line)
                result = self.cursor.fetchone()
                if result != None:
                    continue_check = True
            continue_check = self.comm.mpi_bcast(continue_check)
            if continue_check:
                continue
            red = pf.current_redshift
            # Read the halos off the disk using the Halo Profiler tools.
            hp = HaloProfiler(file, halo_list_file='MergerHalos.out',
                              halo_list_format={'id':0, 'mass':1, 'numpart':2, 'center':[7, 8, 9], 'velocity':[10, 11, 12], 'r_max':13})
            if len(hp.all_halos) == 0:
                mylog.info("Dataset %s has no halos." % file)
                self.with_halos[cycle] = False
                del hp
                continue
            mylog.info("Entering halos into database for z=%f" % red)
            if self.comm.rank == 0:
                for ID,halo in enumerate(hp.all_halos):
                    numpart = int(halo['numpart'])
                    values = (None, currt, red, ID, halo['mass'], numpart,
                    halo['center'][0], halo['center'][1], halo['center'][2],
                    halo['velocity'][0], halo['velocity'][1], halo['velocity'][2],
                    halo['r_max'] / pf['mpc'],
                    -1,0.,-1,0.,-1,0.,-1,0.,-1,0.)
                    # 23 question marks for 23 data columns.
                    line = ''
                    for i in range(23):
                        line += '?,'
                    # Pull off the last comma.
                    line = 'INSERT into Halos VALUES (' + line[:-1] + ')'
                    self.cursor.execute(line, values)
                self.conn.commit()
            self.comm.barrier()
            del hp
    
    def _open_create_database(self):
        # open the database. This creates the database file on disk if it
        # doesn't already exist. Open it on root only.
        self.conn = sql.connect(self.database)
        self.cursor = self.conn.cursor()

    def _create_halo_table(self):
        # Handle the error if the table already exists by doing nothing.
        try:
            self.cursor.execute(create_db_line)
            self.conn.commit()
        except sql.OperationalError:
            pass
    
    def _find_likely_children(self, parentfile, childfile):
        # For each halo in the parent list, identify likely children in the 
        # list of children.

        # First, read in the locations of the child halos.
        child_pf = load(childfile)
        child_t = child_pf.unique_identifier
        if self.comm.rank == 0:
            line = "SELECT SnapHaloID, CenMassX, CenMassY, CenMassZ FROM \
            Halos WHERE SnapCurrentTimeIdentifier = %d" % child_t
            self.cursor.execute(line)
            
            mylog.info("Finding likely parents for z=%1.5f child halos." % \
                child_pf.current_redshift)
            
            # Build the kdtree for the children by looping over the fetched rows.
            # Normalize the points for use only within the kdtree.
            child_points = []
            for row in self.cursor:
                child_points.append([row[1] / self.period[0],
                row[2] / self.period[1],
                row[3] / self.period[2]])
            child_points = na.array(child_points)
            kdtree = cKDTree(child_points, leafsize = 10)
    
        # Find the parent points from the database.
        parent_pf = load(parentfile)
        parent_t = parent_pf.unique_identifier
        if self.comm.rank == 0:
            line = "SELECT SnapHaloID, CenMassX, CenMassY, CenMassZ FROM \
            Halos WHERE SnapCurrentTimeIdentifier = %d" % parent_t
            self.cursor.execute(line)
    
            # Loop over the returned rows, and find the likely neighbors for the
            # parents.
            candidates = {}
            for row in self.cursor:
                # Normalize positions for use within the kdtree.
                query = na.array([row[1] / self.period[0],
                row[2] / self.period[1],
                row[3] / self.period[2]])
                NNtags = kdtree.query(query, NumNeighbors, period=self.period)[1]
                nIDs = []
                for n in NNtags:
                    if n not in nIDs:
                        nIDs.append(n)
                # We need to fill in fake halos if there aren't enough halos,
                # which can happen at high redshifts.
                while len(nIDs) < NumNeighbors:
                    nIDs.append(-1)
                candidates[row[0]] = nIDs
            del kdtree
        else:
            candidates = None

        # Sync across tasks.
        candidates = self.comm.mpi_bcast(candidates)
        self.candidates = candidates
        
        # This stores the masses contributed to each child candidate.
        # The +1 is an extra element in the array that collects garbage
        # values. This is allowing us to eliminate a try/except later.
        # This extra array element will be cut off eventually.
        self.child_mass_arr = na.zeros(len(candidates)*NumNeighbors + 1,
            dtype='float64')
        # Records where to put the entries in the above array.
        self.child_mass_loc = defaultdict(dict)
        # Fill it out with sub-nested default dicts that point to the
        # garbage slot, and then fill it will correct values for (possibly)
        # related parent/child halo pairs.
        for i,halo in enumerate(sorted(candidates)):
            self.child_mass_loc[halo] = defaultdict(minus_one)
            for j, child in enumerate(candidates[halo]):
                self.child_mass_loc[halo][child] = i*NumNeighbors + j

    def _build_h5_refs(self, filename):
        # For this snapshot, add lists of file names that contain the
        # particle info for each halo.
        if not hasattr(self, 'h5files'):
            self.h5files = defaultdict(dict)
        if not hasattr(self, 'names'):
            self.names = defaultdict(set)
        file_pf = load(filename)
        currt = file_pf.unique_identifier
        dir = os.path.dirname(filename)
        h5txt = os.path.join(dir, 'MergerHalos.txt')
        lines = file(h5txt)
        names = set([])
        for i,line in enumerate(lines):
            # Get rid of the carriage returns and turn it into a list.
            line = line.strip().split()
            self.h5files[currt][i] = line[1:]
            names.update(line[1:])
            self.names[currt].update(line[1:])
        lines.close()

    def _compute_child_fraction(self, parentfile, childfile, last):
        # Given a parent and child snapshot, and a list of child candidates,
        # compute what fraction of the parent halo goes to each of the children.
        
        parent_pf = load(parentfile)
        child_pf = load(childfile)
        parent_currt = parent_pf.unique_identifier
        child_currt = child_pf.unique_identifier
        
        mylog.info("Computing fractional contribututions of particles to z=%1.5f halos." % \
            child_pf.current_redshift)
        
        if last == None:
            # First we're going to read in the particles, haloIDs and masses from
            # the parent dataset.
            parent_names = list(self.names[parent_currt])
            parent_names.sort()
            parent_IDs = []
            parent_masses = []
            parent_halos = []
            for i,pname in enumerate(parent_names):
                if i>=self.comm.rank and i%self.comm.size==self.comm.rank:
                    h5fp = h5py.File(pname)
                    for group in h5fp:
                        gID = int(group[4:])
                        thisIDs = h5fp[group]['particle_index'][:]
                        thisMasses = h5fp[group]['ParticleMassMsun'][:]
                        parent_IDs.append(thisIDs)
                        parent_masses.append(thisMasses)
                        parent_halos.append(na.ones(len(thisIDs),
                            dtype='int32') * gID)
                        del thisIDs, thisMasses
                    h5fp.close()
            # Sort the arrays by particle index in ascending order.
            if len(parent_IDs)==0:
                parent_IDs = na.array([], dtype='int64')
                parent_masses = na.array([], dtype='float64')
                parent_halos = na.array([], dtype='int32')
            else:
                parent_IDs = na.concatenate(parent_IDs).astype('int64')
                parent_masses = na.concatenate(parent_masses).astype('float64')
                parent_halos = na.concatenate(parent_halos).astype('int32')
                sort = parent_IDs.argsort()
                parent_IDs = parent_IDs[sort]
                parent_masses = parent_masses[sort]
                parent_halos = parent_halos[sort]
                del sort
        else:
            # We can use old data and save disk reading.
            (parent_IDs, parent_masses, parent_halos) = last
        # Used to communicate un-matched particles.
        parent_send = na.ones(parent_IDs.size, dtype='bool')

        # Now get the child halo data.
        child_names = list(self.names[child_currt])
        child_names.sort()
        child_IDs = []
        child_masses = []
        child_halos = []
        for i,cname in enumerate(child_names):
            if i>=self.comm.rank and i%self.comm.size==self.comm.rank:
                h5fp = h5py.File(cname)
                for group in h5fp:
                    gID = int(group[4:])
                    thisIDs = h5fp[group]['particle_index'][:]
                    thisMasses = h5fp[group]['ParticleMassMsun'][:]
                    child_IDs.append(thisIDs)
                    child_masses.append(thisMasses)
                    child_halos.append(na.ones(len(thisIDs),
                        dtype='int32') * gID)
                    del thisIDs, thisMasses
                h5fp.close()
        # Sort the arrays by particle index in ascending order.
        if len(child_IDs)==0:
            child_IDs = na.array([], dtype='int64')
            child_masses = na.array([], dtype='float64')
            child_halos = na.array([], dtype='int32')
        else:
            child_IDs = na.concatenate(child_IDs).astype('int64')
            child_masses = na.concatenate(child_masses)
            child_halos = na.concatenate(child_halos)
            sort = child_IDs.argsort()
            child_IDs = child_IDs[sort]
            child_masses = child_masses[sort]
            child_halos = child_halos[sort]
            del sort

        child_send = na.ones(child_IDs.size, dtype='bool')
        
        # Match particles in halos.
        self._match(parent_IDs, child_IDs, parent_halos, child_halos,
            parent_masses, parent_send, child_send)

        # Now we send all the un-matched particles to the root task for one more
        # pass. This depends on the assumption that most of the particles do
        # not move very much between data dumps, so that not too many particles
        # will be dumped on the single task.
        parent_IDs_tosend = parent_IDs[parent_send]
        parent_masses_tosend = parent_masses[parent_send]
        parent_halos_tosend = parent_halos[parent_send]
        child_IDs_tosend = child_IDs[child_send]
        child_halos_tosend = child_halos[child_send]
        del parent_send, child_send
        
        parent_IDs_tosend = self.comm.par_combine_object(parent_IDs_tosend,
                datatype="array", op="cat")
        parent_masses_tosend = self.comm.par_combine_object(parent_masses_tosend,
                datatype="array", op="cat")
        parent_halos_tosend = self.comm.par_combine_object(parent_halos_tosend,
                datatype="array", op="cat")
        child_IDs_tosend = self.comm.par_combine_object(child_IDs_tosend,
                datatype="array", op="cat")
        child_halos_tosend = self.comm.par_combine_object(child_halos_tosend,
                datatype="array", op="cat")

        # Resort the received particles.
        Psort = parent_IDs_tosend.argsort()
        parent_IDs_tosend = parent_IDs_tosend[Psort]
        parent_masses_tosend = parent_masses_tosend[Psort]
        parent_halos_tosend = parent_halos_tosend[Psort]
        Csort = child_IDs_tosend.argsort()
        child_IDs_tosend = child_IDs_tosend[Csort]
        child_halos_tosend = child_halos_tosend[Csort]
        del Psort, Csort

        # Now again, but only on the root task.
        if self.comm.rank == 0:
            self._match(parent_IDs_tosend, child_IDs_tosend,
            parent_halos_tosend, child_halos_tosend, parent_masses_tosend)

        # Now we sum up the contributions globally.
        self.child_mass_arr = self.comm.mpi_allreduce(self.child_mass_arr)
        
        # Trim off the garbage collection.
        self.child_mass_arr = self.child_mass_arr[:-1]
        
        if self.comm.rank == 0:
            # Turn these Msol masses into percentages of the parent.
            line = "SELECT HaloMass FROM Halos WHERE SnapCurrentTimeIdentifier=%d \
            ORDER BY SnapHaloID ASC;" % parent_currt
            self.cursor.execute(line)
            mark = 0
            result = self.cursor.fetchone()
            while result:
                mass = result[0]
                self.child_mass_arr[mark:mark+NumNeighbors] /= mass
                mark += NumNeighbors
                result = self.cursor.fetchone()
            
            # Get the global ID for the SnapHaloID=0 from the child, this will
            # be used to prevent unnecessary SQL reads.
            line = "SELECT GlobalHaloID FROM Halos WHERE SnapCurrentTimeIdentifier=%d \
            AND SnapHaloID=0;" % child_currt
            self.cursor.execute(line)
            baseChildID = self.cursor.fetchone()[0]
        else:
            baseChildID = None
        
        # Sync up data on all tasks.
        self.child_mass_arr = self.comm.mpi_bcast(self.child_mass_arr)
        baseChildID = self.comm.mpi_bcast(baseChildID)
        
        # Now we prepare a big list of writes to put in the database.
        for i,parent_halo in enumerate(sorted(self.candidates)):
            child_indexes = []
            child_per = []
            for j,child in enumerate(self.candidates[parent_halo]):
                if child == -1:
                    # Account for fake children.
                    child_indexes.append(-1)
                    child_per.append(0.)
                    continue
                # We need to get the GlobalHaloID for this child.
                child_globalID = baseChildID + child
                child_indexes.append(child_globalID)
                child_per.append(self.child_mass_arr[i*NumNeighbors + j])
            # Sort by percentages, desending.
            child_per, child_indexes = zip(*sorted(zip(child_per, child_indexes), reverse=True))
            values = []
            for pair_count, pair in enumerate(zip(child_indexes, child_per)):
                if pair_count == NumDB: break
                values.extend([int(pair[0]), float(pair[1])])
            #values.extend([parent_currt, parent_halo])
            # This has the child ID, child percent listed NumDB times, followed
            # by the currt and this parent halo ID (SnapHaloID).
            #values = tuple(values)
            self.write_values.append(values)
            self.write_values_dict[parent_currt][parent_halo] = values

        # Clean up.
        del parent_IDs, parent_masses, parent_halos
        del parent_IDs_tosend, parent_masses_tosend
        del parent_halos_tosend, child_IDs_tosend, child_halos_tosend
        gc.collect()
        
        return (child_IDs, child_masses, child_halos)

    def _match(self, parent_IDs, child_IDs, parent_halos, child_halos,
            parent_masses, parent_send = None, child_send = None):
        # Pick out IDs that are in both arrays.
        parent_in_child = na.in1d(parent_IDs, child_IDs, assume_unique = True)
        child_in_parent = na.in1d(child_IDs, parent_IDs, assume_unique = True)
        # Pare down the arrays to just matched particle IDs.
        parent_halos_cut = parent_halos[parent_in_child]
        child_halos_cut = child_halos[child_in_parent]
        parent_masses_cut = parent_masses[parent_in_child]
        # Mark the IDs that have matches so they're not sent later.
        if parent_send is not None:
            parent_send[parent_in_child] = False
            child_send[child_in_parent] = False
        # For matching pairs of particles, add the contribution of the mass.
        # Occasionally, there are matches of particle IDs where the parent
        # and child halos have not been identified as likely relations,
        # and in that case loc will be returned as -1, which is the 'garbage'
        # position in child_mass_arr. This will be trimmed off later.
        for i,pair in enumerate(zip(parent_halos_cut, child_halos_cut)):
            loc = self.child_mass_loc[pair[0]][pair[1]]
            self.child_mass_arr[loc] += parent_masses_cut[i]
        if parent_send is None:
            mylog.info("Clean-up round matched %d of %d parents and %d children." % \
            (parent_in_child.sum(), parent_IDs.size, child_IDs.size))

    def _copy_and_update_db(self):
        """
        Because doing an UPDATE of a SQLite database is really slow, what we'll
        do here is basically read in lines from the database, and then insert
        the parent-child relationships, writing to a new DB.
        """
        # All of this happens only on the root task!
        temp_name = self.database + '-tmp'
        to_write = []
        # Open the temporary database.
        try:
            os.remove(temp_name)
        except OSError:
            pass
        temp_conn = sql.connect(temp_name)
        temp_cursor = temp_conn.cursor()
        line = "CREATE TABLE Halos (GlobalHaloID INTEGER PRIMARY KEY,\
                SnapCurrentTimeIdentifier INTEGER, SnapZ FLOAT, SnapHaloID INTEGER, \
                HaloMass FLOAT,\
                NumPart INTEGER, CenMassX FLOAT, CenMassY FLOAT,\
                CenMassZ FLOAT, BulkVelX FLOAT, BulkVelY FLOAT, BulkVelZ FLOAT,\
                MaxRad FLOAT,\
                ChildHaloID0 INTEGER, ChildHaloFrac0 FLOAT, \
                ChildHaloID1 INTEGER, ChildHaloFrac1 FLOAT, \
                ChildHaloID2 INTEGER, ChildHaloFrac2 FLOAT, \
                ChildHaloID3 INTEGER, ChildHaloFrac3 FLOAT, \
                ChildHaloID4 INTEGER, ChildHaloFrac4 FLOAT);"
        temp_cursor.execute(line)
        temp_conn.commit()
        # Get all the data!
        self.cursor.execute("SELECT * FROM Halos;")
        results = self.cursor.fetchone()
        while results:
            results = list(results)
            currt = results[1]
            hid = results[3]
            # If for some reason this halo doesn't have relationships,
            # we'll just keep the old results the same.
            try:
                lookup = self.write_values_dict[currt][hid]
                new = tuple(results[:-10] + lookup)
            except KeyError:
                new = tuple(results)
            to_write.append(new)
            results = self.cursor.fetchone()
        # Now write to the temp database.
        # 23 question marks for 23 data columns.
        line = ''
        for i in range(23):
            line += '?,'
        # Pull off the last comma.
        line = 'INSERT into Halos VALUES (' + line[:-1] + ')'
        for insert in to_write:
            temp_cursor.execute(line, insert)
        temp_conn.commit()
        temp_cursor.close()
        temp_conn.close()
        self._close_database()
        os.rename(temp_name, self.database)

    def _write_index(self):
        mylog.info("Creating database index.")
        line = "CREATE INDEX IF NOT EXISTS HalosIndex ON Halos ("
        for name in columns:
            line += name +","
        line = line[:-1] + ");"
        self.cursor.execute(line)

    def _find_progress(self):
        # This queries the database to see how far along work has already come
        # to identify parent->child relationships.
        line = """SELECT ChildHaloID0, SnapZ from halos WHERE SnapHaloID = 0
        ORDER BY SnapZ DESC;"""
        self.cursor.execute(line)
        results = self.cursor.fetchone()
        while results:
            results = list(results)
            if results[0] == -1:
                # We've hit a dump that does not have relationships. Save this.
                return results[1] # the SnapZ.
            results = self.cursor.fetchone()
        return 0.

class MergerTreeConnect(DatabaseFunctions):
    def __init__(self, database='halos.db'):
        r"""Create a convenience object for accessing data from the halo database.
        
        Parameters
        ----------
        database : String
            The name of the halo database to access. Default = 'halos.db'.
        
        Examples
        -------
        >>> mtc = MergerTreeConnect('/home/user/sim1-halos.db')
        """
        self.database = database
        result = self._open_database()
        if not result:
            return None
    
    def close(self):
        r"""Cleanly close access to the database.
        
        Examples
        --------
        >>> mtc.close()
        """
        # To be more like typical Python open/close.
        self._close_database()
    
    def query(self, string):
        r"""Performs a query of the database and returns the results as a list
        of tuple(s), even if the result is singular.
        
        Parameters
        ----------
        string : String
            The SQL query of the database.
        
        Examples
        -------
        >>> results = mtc.query("SELECT GlobalHaloID from Halos where SnapHaloID = 0 and \
        ... SnapZ = 0;")
        """
        # Query the database and return a list of tuples.
        if string is None:
            mylog.error("You must enter a SQL query.")
            return None
        items = []
        self.cursor.execute(string)
        results = self.cursor.fetchone()
        while results:
            items.append(results)
            results = self.cursor.fetchone()
        return items

    def get_GlobalHaloID(self, SnapHaloID, z):
        r"""Returns the GlobalHaloID for the given halo.
        
        Parameters
        ---------
        SnapHaloID : Integer
            The index label for the halo of interest, equivalent to
            the first column of the halo finder text output file.
        z : Float
            The redshift for the halo of interest. The value returned will be
            for the halo with SnapHaloID equal to ID (above) with redshift
            closest to this value.
        
        Examples
        --------
        >>> this_halo = mtc.get_GlobalHaloID(0, 0.)
        """
        string = "SELECT GlobalHaloID,SnapZ FROM Halos WHERE SnapHaloID = %d;" \
            % SnapHaloID
        minz = 99999.
        # If -1 is returned, something went wrong.
        this_halo = -1
        self.cursor.execute(string)
        results = self.cursor.fetchone()
        while results:
            if abs(results[1] - z) < minz:
                minz = abs(results[1] - z)
                this_halo = results[0]
            results = self.cursor.fetchone()
        return this_halo

    def get_halo_parents(self, GlobalHaloID):
        r"""Returns a list of the parent halos to the given halo, along with
        the contribution fractions from parent to child.
        
        This function returns a list of lists, where each entry in the top list
        is [GlobalHaloID, ChildHaloFrac] of the parent halo in relationship
        to the given child halo.
        
        Parameters
        ----------
        GlobalHaloID : Integer
            The GlobalHaloID of the halo of interest.
        
        Examples
        --------
        >>> parents = mtc.get_halo_parents(1688)
        >>> print parents
        [[1544, 0.9642857141249418],
         [1613, 0.0],
         [1614, 0.0],
         [1489, 0.0],
         [1512, 0.0],
         [1519, 0.0],
         [1609, 0.0]]
        """
        parents = []
        for i in range(NumDB):
            string = "SELECT GlobalHaloID, ChildHaloFrac%d FROM Halos\
            WHERE ChildHaloID%d=%d;" % (i, i, GlobalHaloID)
            self.cursor.execute(string)
            results = self.cursor.fetchone()
            while results:
                parents.append([results[0], results[1]])
                results = self.cursor.fetchone()
        return parents

    def get_direct_parent(self, GlobalHaloID):
        r"""Returns the GlobalHaloID of the direct parent of the given halo.
        
        This is accomplished by identifying the most massive parent halo
        that contributes at least 50% of its mass to the given halo.
        
        Parameters
        ----------
        GlobalHaloID : Integer
            The GlobalHaloID of the halo of interest.
        
        Examples
        --------
        >>> parent = mtc.get_direct_parent(1688)
        >>> print parent
        1544
        """
        parents = self.get_halo_parents(GlobalHaloID)
        mass = 0
        ID = None
        for parent in parents:
            if parent[1] < 0.5: continue
            info = self.get_halo_info(parent[0])
            if info['HaloMass'] > mass:
                mass = info['HaloMass']
                ID = parent[0]
        return ID

    def get_halo_info(self, GlobalHaloID):
        r"""Returns all available information for the given GlobalHaloID
        in the form of a dict.
        
        Parameters
        ----------
        GlobalHaloID : Integer
            The unique index for the halo of interest.
        
        Examples
        --------
        >>> info = mtc.get_halo_info(1544)
        >>> print info
        {'BulkVelX': -32759799.359999999,
         'BulkVelY': -28740239.109999999,
         'BulkVelZ': -20066000.690000001,
         'CenMassX': 0.23059111360000001,
         'CenMassY': 0.4061139809,
         'CenMassZ': 0.80882763749999997,
         'ChildHaloFrac0': 0.9642857141249418,
         'ChildHaloFrac1': 0.0,
         'ChildHaloFrac2': 0.0,
         'ChildHaloFrac3': 0.0,
         'ChildHaloFrac4': 0.0,
         'ChildHaloID0': 1688,
         'ChildHaloID1': 1712,
         'ChildHaloID2': 1664,
         'ChildHaloID3': 1657,
         'ChildHaloID4': 1634,
         'GlobalHaloID': 1544,
         'HaloMass': 20934692770000.0,
         'MaxRad': 0.01531299899,
         'NumPart': 196,
         'SnapCurrentTimeIdentifier': 1275946788,
         'SnapHaloID': 56,
         'SnapZ': 0.024169713061444002}
        """
        string = "SELECT * FROM Halos WHERE GlobalHaloID=%d;" % GlobalHaloID
        d = {}
        self.cursor.execute(string)
        results = self.cursor.fetchone()
        for pair in zip(columns, results):
            d[pair[0]] = pair[1]
        return d

class Node(object):
    def __init__(self, CoM, mass, parentIDs, z, color):
        self.CoM = CoM
        self.mass = mass
        self.parentIDs = parentIDs # In descending order of contribution
        self.z = z
        self.color = color

class Link(object):
    def __init__(self):
        self.childIDs = []
        self.fractions = []

class MergerTreeDotOutput(DatabaseFunctions, ParallelAnalysisInterface):
    def __init__(self, halos=None, database='halos.db',
            dotfile='MergerTree.gv', current_time=None, link_min=0.2):
        r"""Output the merger tree history for a given set of halo(s) in Graphviz
        format.
        
        Parameters
        ---------
        halos : Integer or list of integers
            If current_time below is not specified or is None, this is an integer
            or list of integers with the GlobalHaloIDs of the halos to be
            tracked. If current_time is specified, this is the SnapHaloIDs
            for the halos to be tracked, which is identical to what is in
            HopAnalysis.out files (for example).
        database : String
            The name of the database file. Default = 'halos.db'.
        dotfile : String
            The name of the file to write to. Default = 'MergerTree.gv'.
            The suffix of this name gives the format of the output file,
            so 'MergerTree.jpg' would output a jpg file. "dot -v" (from the
            command line) will print
            a list of image formats supported on the system. The default
            suffix '.gv' will output the results to a text file in the Graphviz
            markup language.
        current_time : Integer
            The SnapCurrentTimeIdentifier for the snapshot for the halos to
            be tracked. This is identical to the CurrentTimeIdentifier in
            Enzo restart files. Default = None.
        link_min : Float
            When establishing a parent/child relationship, this is the minimum
            mass fraction of the parent halo contributed to
            the child halo that will be tracked
            while building the Graphviz file. Default = 0.2.
        
        Examples
        --------
        >>> MergerTreeDotOutput(halos=182842, database='/home/user/sim1-halos.db',
        ... dotfile = 'halo-182842.gv')
        """
        ParallelAnalysisInterface.__init__(self)
        self.database = database
        self.link_min = link_min
        if halos is None:
            mylog.error("Please provide at least one halo to start the tree. Exiting.")
            return None
        result = self._open_database()
        if not result:
            mylog.warn("The database did not open correctly!")
            return None
        if type(halos) == types.IntType:
            halos = [halos]
        if current_time is not None:
            halos = self._translate_haloIDs(halos, current_time)
        newhalos = set(halos)
        # Create the pydot graph object.
        self.graph = pydot.Dot('galaxy', graph_type='digraph')
        # Build some initially empty subgraphs, which are used to identify
        # nodes that are on the same rank (redshift).
        line = "SELECT DISTINCT SnapZ FROM Halos;"
        self.cursor.execute(line)
        self.subgs = {}
        result = self.cursor.fetchone()
        while result:
            self.subgs[result[0]] = pydot.Subgraph('', rank = 'same')
            self.graph.add_subgraph(self.subgs[result[0]])
            result = self.cursor.fetchone()
        # For the first set of halos.
        self._add_nodes(newhalos)
        # Recurse over parents.
        while len(newhalos) > 0:
            mylog.info("Finding parents for %d children." % len(newhalos))
            newhalos = self._find_parents(newhalos)
            self._add_nodes(newhalos)
        self._write_dotfile(dotfile)
        return None

    def _translate_haloIDs(self, halos, current_time):
        # If the input is in the haloID equivalent to SnapHaloID, translate them
        # to GlobalHaloIDs.
        new_haloIDs=[]
        for halo in halos:
            line = "SELECT GlobalHaloID FROM Halos WHERE SnapHaloID=? AND \
            SnapCurrentTimeIdentifier=? limit 1;"
            values = (halo, current_time)
            self.cursor.execute(line, values)
            new_haloIDs.append(self.cursor.fetchone()[0])
        return new_haloIDs
        
    def _find_parents(self, halos):
        # Given a set of halos, find their parents and add that to each of their
        # node records. At the same time, make a link record for that
        # relationship.
        # This stores the newly discovered parent halos.
        newhalos = set([])
        for halo in halos:
            line = "SELECT GlobalHaloID, ChildHaloFrac0,\
                ChildHaloFrac1, ChildHaloFrac2,ChildHaloFrac3, ChildHaloFrac4,\
                ChildHaloID0, ChildHaloID1, ChildHaloID2, \
                ChildHaloID3, ChildHaloID4 \
                FROM Halos WHERE\
                ChildHaloID0=? or ChildHaloID1=? or ChildHaloID2=? or\
                ChildHaloID3=? or ChildHaloID4=?;"
            values = (halo, halo, halo, halo, halo)
            self.cursor.execute(line, values)
            result = self.cursor.fetchone()
            while result:
                res = list(result)
                pID = result[0]
                pfracs = res[1:6]
                cIDs = res[6:11]
                for pair in zip(cIDs, pfracs):
                    if pair[1] <= self.link_min or pair[0] != halo:
                        continue
                    else:
                        self.graph.add_edge(pydot.Edge(pID, halo,
                        label = "%3.2f%%" % float(pair[1]*100),
                        color = "blue", 
                        fontsize = "10"))
                        newhalos.add(pID)
                result = self.cursor.fetchone()
        return newhalos
    
    def _add_nodes(self, newhalos):
        # Each call of this function always happens for a set of newhalos that
        # are at the same z. To give the halos color we will figure out how
        # many halos total were found this z.
        # There's probably a way to do this with only one SQL operation.
        if len(newhalos) == 0:
            return
        ahalo = list(newhalos)[0]
        line = 'SELECT SnapCurrentTimeIdentifier FROM Halos WHERE GlobalHaloID=?;'
        values = (ahalo,)
        self.cursor.execute(line, values)
        result = self.cursor.fetchone()
        # Use currt to get the number.
        line = 'SELECT max(SnapHaloID) FROM Halos where SnapCurrentTimeIdentifier=?;'
        values = (result[0],)
        self.cursor.execute(line, values)
        maxID = self.cursor.fetchone()[0]
        # For the new halos, create nodes for them.
        for halo in newhalos:
            line = 'SELECT SnapZ, HaloMass, CenMassX, CenMassY, CenMassZ,\
            SnapHaloID FROM Halos WHERE GlobalHaloID=? limit 1;'
            value = (halo,)
            self.cursor.execute(line, value)
            result = self.cursor.fetchone()
            # Add the node to the pydot graph.
            color_float = 1. - float(result[5])/(maxID+1)
            self.graph.add_node(pydot.Node(halo,
                label = "{%1.3e\\n(%1.3f,%1.3f,%1.3f)}" % \
                (result[1], result[2], result[3], result[4]),
                shape = "record",
                color = "%0.3f 1. %0.3f" % (color_float, color_float)))
            # Add this node to the correct subgraph.
            self.subgs[result[0]].add_node(pydot.Node(halo))
            # If this was the first node added to this subgraph, also add
            # the lone node for the redshift value.
            if len(self.subgs[result[0]].get_node_list()) == 1:
                self.subgs[result[0]].add_node(pydot.Node("%1.5e" % result[0],
                label = "%1.5f" % result[0],
                shape = "record", color = "green"))

    def _write_dotfile(self, dotfile):
        # Based on the suffix of the file name, write out the result to a file.
        suffix = dotfile.split(".")[-1]
        if suffix == "gv": suffix = "raw"
        mylog.info("Writing %s format %s to disk." % (suffix, dotfile))
        self.graph.write("%s" % dotfile, format=suffix)

class MergerTreeTextOutput(DatabaseFunctions, ParallelAnalysisInterface):
    def __init__(self, database='halos.db', outfile='MergerTreeDB.txt'):
        r"""Dump the contents of the merger tree database to a text file.
        This is generally not recommended.
        
        Parameters
        ----------
        database : String
            Name of the database to access. Default = 'halos.db'.
        outfile : String
            Name of the file to write to. Default = 'MergerTreeDB.txt'.
        
        Examples
        --------
        >>> MergerTreeTextOutput(database='/home/user/sim1-halos.db',
        ... outfile='halos-db.txt')
        """
        ParallelAnalysisInterface.__init__(self)
        self.database = database
        self.outfile = outfile
        result = self._open_database()
        if not result:
            mylog.warn("Database file not read correctly!")
            return None
        self._write_out()
        self._close_database()
        return None
    
    def _write_out(self):
        # Essentially dump the contents of the database into a text file.
        fp = open(self.outfile, "w")
        # Make the header line.
        spacing = {}
        for column in columns:
            spacing[column] = (max(15,len(column)+1))
        line = "# "
        for column in columns:
            line += "%s" % column.ljust(spacing[column])
        line += "\n"
        fp.write(line)
        # Get the data.
        line = "SELECT * FROM Halos ORDER BY SnapZ DESC, SnapHaloID ASC;"
        self.cursor.execute(line)
        results = self.cursor.fetchone()
        # Write out the columns.
        while results:
            line = "  "
            for i,column in enumerate(columns):
                if column_types[column] == "FLOAT":
                    this = "%1.6e" % results[i]
                    line += this.ljust(spacing[column])
                if column_types[column] == "INTEGER":
                    this = "%d" % results[i]
                    line += this.ljust(spacing[column])
            line += "\n"
            fp.write(line)
            results = self.cursor.fetchone()
        fp.close()
        

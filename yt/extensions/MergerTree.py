"""
MergerTree class and member functions.

Author: Stephen Skory <sskory@physics.ucsd.edu>
Affiliation: CASS/UC San Diego, CA
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2010 Stephen Skory.  All Rights Reserved.

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

import yt.lagos as lagos
from yt.lagos.kd import *
from yt.lagos.HaloFinding import HaloFinder
from yt.logger import lagosLogger as mylog
import yt.extensions.HaloProfiler as HP

import numpy as na
import os, glob
import h5py
import types
import sqlite3 as sql
from collections import defaultdict

class MergerTree(lagos.ParallelAnalysisInterface):
    def __init__(self, restart_files=[], database='halos.db',
            halo_finder_function=HaloFinder, halo_finder_threshold=80.0):
        self.restart_files = restart_files # list of enzo restart files
        self.database = database # the sqlite database of haloes.
        self.halo_finder_function = halo_finder_function # which halo finder to use
        self.halo_finder_threshold = halo_finder_threshold # overdensity threshold
        self._open_database()
        self._create_halo_table()
        self._run_halo_finder_add_to_db()
        # Find the h5 file names for all the halos.
        for snap in self.restart_files:
            self._build_h5_refs(snap)
        # Loop over the pairs of snapshots to locate likely neighbors, and
        # then use those likely neighbors to compute fractional contributions.
        for pair in zip(self.restart_files[:-1], self.restart_files[1:]):
            candidates = self._find_likely_children(pair[0], pair[1])
            self._compute_child_fraction(pair[0], pair[1], candidates)
        self._close_database()
        self._close_h5fp()
        
    def _read_halo_lists(self):
        self.halo_lists = []
        for i,file in enumerate(self.halo_files):
            hp = HP.HaloProfiler(self.restart_files[i], halo_list_file=file)
            self.halo_lists.append(hp.all_halos)

    def _run_halo_finder_add_to_db(self):
        for file in self.restart_files:
            pf = lagos.EnzoStaticOutput(file)
            # If the halos are already found, skip this one.
            dir = os.path.dirname(file)
            if os.path.exists(os.path.join(dir, 'MergerHalos.out')) and \
                    os.path.exists(os.path.join(dir, 'MergerHalos.txt')) and \
                    glob.glob(os.path.join(dir, 'MergerHalos*h5')) is not []:
                pass
            else:
                # Run the halo finder.
                halos = self.halo_finder_function(pf,
                    threshold=self.halo_finder_threshold, dm_only=True)
                halos.write_out(os.path.join(dir, 'MergerHalos.out'))
                halos.write_particle_lists(os.path.join(dir, 'MergerHalos'))
                halos.write_particle_lists_txt(os.path.join(dir, 'MergerHalos'))
                del halos
            # Now add halo data to the db if it isn't already there by
            # checking the first halo.
            currt = pf['CurrentTimeIdentifier']
            line = "SELECT GlobalHaloID from Halos where SnapHaloID=0\
            and SnapCurrentTimeIdentifier=%d;" % currt
            self.cursor.execute(line)
            result = self.cursor.fetchone()
            if result != None:
                continue
            red = pf['CosmologyCurrentRedshift']
            # Read the halos off the disk using the Halo Profiler tools.
            hp = HP.HaloProfiler(file, halo_list_file='MergerHalos.out',
            halo_list_format={'id':0, 'mass':1, 'numpart':2, 'center':[7, 8, 9], 'velocity':[10, 11, 12], 'r_max':13})
            mylog.info("Entering halos into database for z=%f" % red)
            for ID,halo in enumerate(hp.all_halos):
                numpart = int(halo['numpart'])
                values = (None, currt, red, ID, halo['mass'], numpart,
                halo['center'][0], halo['center'][1], halo['center'][2],
                halo['velocity'][0], halo['velocity'][1], halo['velocity'][2],
                halo['r_max'],
                -1,0.,-1,0.,-1,0.,-1,0.,-1,0.)
                # 23 question marks for 23 data columns.
                line = ''
                for i in range(23):
                    line += '?,'
                # Pull off the last comma.
                line = 'INSERT into Halos VALUES (' + line[:-1] + ')'
                self.cursor.execute(line, values)
            self.conn.commit()
            del hp
    
    def _open_database(self):
        # open the database. This creates the database if it doesn't already exist.
        self.conn = sql.connect(self.database)
        self.cursor = self.conn.cursor()
    
    def _close_database(self):
        # close the database cleanly.
        self.conn.commit()
        self.cursor.close()

    def _create_halo_table(self):
        # Handle the error if it already exists.
        try:
            # Create the table that will store the halo data.
            line = "CREATE TABLE Halos (GlobalHaloID INTEGER PRIMARY KEY,\
                SnapCurrentTimeIdentifier INTEGER, SnapZ FLOAT, SnapHaloID INTEGER, \
                DarkMatterMass FLOAT,\
                NumPart INTEGER, CenMassX FLOAT, CenMassY FLOAT,\
                CenMassZ FLOAT, BulkVelX FLOAT, BulkVelY FLOAT, BulkVelZ FLOAT,\
                MaxRad FLOAT,\
                ChildHaloID0 INTEGER, ChildHaloFrac0 FLOAT, \
                ChildHaloID1 INTEGER, ChildHaloFrac1 FLOAT, \
                ChildHaloID2 INTEGER, ChildHaloFrac2 FLOAT, \
                ChildHaloID3 INTEGER, ChildHaloFrac3 FLOAT, \
                ChildHaloID4 INTEGER, ChildHaloFrac4 FLOAT);"
            self.cursor.execute(line)
        except sql.OperationalError:
            return
    
    def _find_likely_children(self, parentfile, childfile):
        # For each halo in the parent list, identify likely children in the 
        # list of children.
        
        # First, read in the locations of the child halos.
        child_pf = lagos.EnzoStaticOutput(childfile)
        child_t = child_pf['CurrentTimeIdentifier']
        line = "SELECT SnapHaloID, CenMassX, CenMassY, CenMassZ FROM \
        Halos WHERE SnapCurrentTimeIdentifier = %d" % child_t
        self.cursor.execute(line)
        
        # Build the kdtree for the children by looping over the fetched rows.
        child_points = []
        for row in self.cursor:
            p = Point()
            p.data = [row[1], row[2], row[3]]
            p.ID = row[0]
            child_points.append(p)
        child_kd = buildKdHyperRectTree(child_points[:],10)

        # Make these just once.
        neighbors = Neighbors()
        neighbors.k = 5

        # Find the parent points from the database.
        parent_pf = lagos.EnzoStaticOutput(parentfile)
        parent_t = parent_pf['CurrentTimeIdentifier']
        line = "SELECT SnapHaloID, CenMassX, CenMassY, CenMassZ FROM \
        Halos WHERE SnapCurrentTimeIdentifier = %d" % parent_t
        self.cursor.execute(line)

        # Loop over the returned rows, and find the likely neighbors for the
        # parents.
        candidates = {}
        for row in self.cursor:
            neighbors.points = []
            neighbors.minDistanceSquared = 100. # should make this a function of the simulation
            cm = [row[1], row[2], row[3]]
            getKNN(cm, child_kd, neighbors, 0., [1.]*3)
            nIDs = []
            for n in neighbors.points:
                nIDs.append(n[1].ID)
            candidates[row[0]] = nIDs
        
        return candidates

    def _build_h5_refs(self, filename):
        # For this snapshot, add lists of file names that contain the
        # particle info for each halo.
        if not hasattr(self, 'h5files'):
            self.h5files = defaultdict(dict)
        file_pf = lagos.EnzoStaticOutput(filename)
        currt = file_pf['CurrentTimeIdentifier']
        dir = os.path.dirname(filename)
        h5txt = os.path.join(dir, 'MergerHalos.txt')
        lines = file(h5txt)
        names = set([])
        for i,line in enumerate(lines):
            # Get rid of the carriage returns and turn it into a list.
            line = line.strip().split()
            self.h5files[currt][i] = line[1:]
            names.update(line[1:])
        lines.close()
        # As long as we're here, open each of the h5 files and store the
        # pointer in a persistent dict.
        if not hasattr(self, 'h5fp'):
            self.h5fp = defaultdict(dict)
        for name in names:
            self.h5fp[currt][name] = h5py.File(name)

    def _close_h5fp(self):
        # Cleanly close the open h5 file pointers.
        for currt in self.h5fp:
            for name in self.h5fp[currt]:
                self.h5fp[currt][name].close()
        
    def _compute_child_fraction(self, parentfile, childfile, candidates):
        # Given a parent and child snapshot, and a list of child candidates,
        # compute what fraction of the parent halo goes to each of the children.
        
        parent_pf = lagos.EnzoStaticOutput(parentfile)
        child_pf = lagos.EnzoStaticOutput(childfile)
        parent_currt = parent_pf['CurrentTimeIdentifier']
        child_currt = child_pf['CurrentTimeIdentifier']
        
        child_percents = {}
        for halo in candidates:
            # Read in its particle IDs
            parent_IDs = na.array([], dtype='int64')
            for h5name in self.h5files[parent_currt][halo]:
                # Get the correct time dict entry, and then the correct h5 file
                # from that snapshot, and then choose this parent's halo
                # group, and then the particle IDs. How's that for a long reach?
                new_IDs = self.h5fp[parent_currt][h5name]['Halo%08d' % halo]['particle_index']
                parent_IDs = na.concatenate((parent_IDs, new_IDs[:]))
            # Loop over its children.
            temp_percents = {}
            for child in candidates[halo]:
                child_IDs = na.array([], dtype='int64')
                for h5name in self.h5files[child_currt][child]:
                    new_IDs = self.h5fp[child_currt][h5name]['Halo%08d' % child]['particle_index']
                    child_IDs = na.concatenate((child_IDs, new_IDs[:]))
                # The IDs shared by both halos.
                intersect = na.intersect1d(parent_IDs, child_IDs)
                # The fraction shared from the parent halo.
                temp_percents[child] = float(intersect.size) / parent_IDs.size
            child_percents[halo] = temp_percents
        
        # Now we write these percents to the existing entry already in the database.
        for halo in candidates:
            child_IDs = []
            child_per = []
            for child in candidates[halo]:
                # We need to get the GlobalHaloID for this child.
                line = 'SELECT GlobalHaloID FROM Halos WHERE \
                SnapCurrentTimeIdentifier=? AND SnapHaloID=?'
                values = (child_currt, child)
                self.cursor.execute(line, values)
                child_globalID = self.cursor.fetchone()[0]
                child_IDs.append(child_globalID)
                child_per.append(child_percents[halo][child])
            # Sort by percentages, desending.
            child_per, child_IDs = zip(*sorted(zip(child_per, child_IDs), reverse=True))
            values = []
            for pair in zip(child_IDs, child_per):
                values.extend(pair)
            values.extend([parent_currt, halo])
            # This has the child ID, child percent listed five times, followed
            # by the currt and this parent halo ID (SnapHaloID).
            values = tuple(values)
            line = 'UPDATE Halos SET ChildHaloID0=?, ChildHaloFrac0=?,\
            ChildHaloID1=?, ChildHaloFrac1=?,\
            ChildHaloID2=?, ChildHaloFrac2=?,\
            ChildHaloID3=?, ChildHaloFrac3=?,\
            ChildHaloID4=?, ChildHaloFrac4=?\
            WHERE SnapCurrentTimeIdentifier=? AND SnapHaloID=?;'
            self.cursor.execute(line, values)
        self.conn.commit()

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

class MergerTreeDotOutput(object):
    def __init__(self, halos, database='halos.db',
            dotfile='MergerTree.dot'):
        if type(halos) == types.IntType:
            halos = [halos]
        newhalos = set(halos)
        self.database = database
        if len(newhalos) == 0:
            mylog.error("Please provide at least one halo to start the tree. Exiting.")
            return None
        result = self._open_database()
        if not result:
            return None
        # A key is the GlobalHaloID for this halo, and the content is a
        # Node object.
        self.nodes = {}
        # A key is the GlobalHaloID for the parent in the relationship,
        # and the content is a Link ojbect.
        self.links = defaultdict(Link)
        # Record which halos are at the same z level for convenience.
        # They key is a z value, and the content a list of co-leveled halo IDs.
        self.levels = defaultdict(list)
        # For the first set of halos.
        self._add_nodes(newhalos)
        # Recurse over parents.
        while len(newhalos) > 0:
            newhalos = self._find_parents(newhalos)
            self._add_nodes(newhalos)
        self._open_dot(dotfile)
        self._write_nodes()
        self._write_links()
        self._write_levels()
        self._close_dot()
        self._close_database()

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
        self.cursor.close()
        
    def _find_parents(self, halos):
        # Given a set of halos, find their parents and add that to each of their
        # node records. At the same time, make a link record for that
        # relationship.
        # This stores the newly discovered parent halos.
        newhalos = set([])
        for halo in halos:
            line = 'SELECT GlobalHaloID, ChildHaloFrac0 from Halos\
            where ChildHaloID0=? or ChildHaloID1=? or ChildHaloID2=?\
            or ChildHaloID3=? or ChildHaloID4=?'
            values = (halo,halo,halo,halo,halo)
            self.cursor.execute(line, values)
            result = self.cursor.fetchone()
            while result:
                pID = result[0]
                pfrac = result[1]
                # Store this.
                self.nodes[halo].parentIDs.append(pID)
                self.links[pID].childIDs.append(halo)
                self.links[pID].fractions.append(pfrac)
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
            line = 'SELECT SnapZ, DarkMatterMass, CenMassX, CenMassY, CenMassZ,\
            SnapHaloID FROM Halos WHERE GlobalHaloID=? limit 1;'
            value = (halo,)
            self.cursor.execute(line, value)
            result = self.cursor.fetchone()
            self.nodes[halo] = Node(na.array([result[2],result[3],result[4]]),
                result[1], [], result[0], 1. - float(result[5])/maxID)
            self.levels[result[0]].append(halo)

    def _open_dot(self, dotfile):
        # Write out the opening stuff in the dotfile.
        self.dotfile=open(dotfile, 'w')
        line = 'digraph galaxy {size="10, 10";\n'
        line += 'node [style=bold, shape=record];\n'
        self.dotfile.write(line)
    
    def _close_dot(self):
        self.dotfile.write("\n};\n")
        self.dotfile.close()
    
    def _write_nodes(self):
        # Write out the nodes to the dot file.
        self.dotfile.write("{\n")
        for halo in self.nodes:
            this = self.nodes[halo]
            line = '"%d" [label="{%1.3e\\n(%1.3f,%1.3f,%1.3f)}", shape="record",' \
                % (halo, this.mass, this.CoM[0], this.CoM[1], this.CoM[2])
            line += ' color="%0.3f 1. %0.3f"];\n' % (this.color, this.color)
            self.dotfile.write(line)
        self.dotfile.write("};\n")
    
    def _write_links(self):
        # Write out the links to the dot file.
        self.dotfile.write("{\n")
        for parent in self.links:
            this = self.links[parent]
            for child,frac in zip(this.childIDs, this.fractions):
                if frac > 0.2:
                    line = '"%d"->"%d" [label="%3.2f%%", color="blue", fontsize=10];\n' \
                        % (parent, child, frac*100.)
                    self.dotfile.write(line)
        self.dotfile.write("};\n")

    def _write_levels(self):
        # Write out the co-leveled halos to the dot file.
        for z in self.levels:
            this = self.levels[z]
            self.dotfile.write("{ rank = same;\n")
            line = '"%1.5f"; ' % z
            for halo in this:
                line += '"%d"; ' % halo
            line += "\n};\n"
            self.dotfile.write(line)
        # Also write out the unlinked boxes for the redshifts.
        line = '{"%1.5f" [label="{%1.5f}", shape="record" color="green"];}\n' \
            % (z, z)
        
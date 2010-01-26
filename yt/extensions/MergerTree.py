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
        self._close_database()
#         self.halo_lists[0] = self._find_likely_children(self.halo_lists[0], self.halo_lists[1])
#         self.halo_lists[0] = self._read_particle_ids(self.halo_lists[0], h5_txts[0])
#         self.halo_lists[1] = self._read_particle_ids(self.halo_lists[1], h5_txts[1])
#         self._find_child_fraction(self.halo_lists[0], self.halo_lists[1])
        
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
                halos = self.halo_finder_function(pf,
                    threshold=self.halo_finder_threshold, dm_only=True)
                halos.write_out(os.path.join(dir, 'MergerHalos.out'))
                halos.write_particle_lists(os.path.join(dir, 'MergerHalos'))
                halos.write_particle_lists_txt(os.path.join(dir, 'MergerHalos'))
                del halos
            # Now add halo data to the db if it isn't already there by
            # checking the first halo.
            currt = pf['CurrentTimeIdentifier']
            line = "SELECT GlobalHaloID from HALOS where SnapHaloID=0\
            and SnapCurrentTimeIdentifier=%d;" % currt
            self.cursor.execute(line)
            result = self.cursor.fetchone()
            if result != None:
                continue
            red = pf['CosmologyCurrentRedshift']
            hp = HP.HaloProfiler(file, halo_list_file='MergerHalos.out',
            halo_list_format={'id':0, 'mass':1, 'numpart':2, 'center':[7, 8, 9], 'velocity':[10, 11, 12], 'r_max':13})
            for ID,halo in enumerate(hp.all_halos):
                numpart = int(halo['numpart'])
                values = (None, currt, red, ID, halo['mass'], numpart,
                halo['center'][0], halo['center'][1], halo['center'][2],
                halo['velocity'][0], halo['velocity'][1], halo['velocity'][2],
                halo['r_max'])
                self.cursor.execute('INSERT into HALOS VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)', values)
            self.conn.commit()
    
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
            line = "CREATE TABLE HALOS (GlobalHaloID INTEGER PRIMARY KEY,\
                SnapCurrentTimeIdentifier INTEGER, SnapZ FLOAT, SnapHaloID INTEGER, \
                DarkMatterMass FLOAT,\
                NumPart INTEGER, CenMassX FLOAT, CenMassY FLOAT,\
                CenMassZ FLOAT, BulkVelX FLOAT, BulkVelY FLOAT, BulkVelZ FLOAT,\
                MaxRad FLOAT);"
            self.cursor.execute(line)
        except sql.OperationalError:
            return
    
    def _find_likely_children(self, parents, children):
        # For each halo in the parent list, identify likely children.
        
        # Build the kdtree for the children.
        child_points = []
        for i,halo in enumerate(children):
            p = Point()
            p.data = halo['center']
            p.ID = i
            child_points.append(p)
        child_kd = buildKdHyperRectTree(child_points[:],10)

        # Make these just once.
        neighbors = Neighbors()
        neighbors.k = self.neighbors

        # Find the closest child halos to this parent.
        for i,parent in enumerate(parents):
            neighbors.points = []
            neighbors.minDistanceSquared = 100. # should make this a function of the simulation
            getKNN(parent['center'], child_kd, neighbors, 0., [1.]*3)
            nIDs = []
            for n in neighbors.points:
                nIDs.append(n[1].ID)
            # Save this to the halo entry.
            parents[i]['neighbors'] = nIDs
        
        return parents
    
    def _read_particle_ids(self, halos, txt_file):
        # Given a txt_file, read in the particles for each halo.
        lines = file(txt_file)
        names = set([])
        for i,line in enumerate(lines):
            # get rid of the carriage returns and turn it into a list
            line = line.strip().split()
            halos[i]['fnames'] = line[1:]
            names.update(line[1:])
        # Open the unique files only once.
        fp_names = {}
        for name in names:
            fp_names[name] = h5py.File(os.path.join(os.path.dirname(txt_file), name))
        # Now read in the indices for all the halos
        for i,halo in enumerate(halos):
            IDs = set([])
            for fname in halo['fnames']:
                tempids = fp_names[fname]['Halo%08d' % i]['particle_index']
                IDs.update(tempids[:])
            # Save the particle IDs
            halos[i]['IDs'] = IDs
        # Close all the various files.
        lines.close()
        for name in names:
            fp_names[name].close()
        return halos
    
    def _find_child_fraction(self, parents, children):
        # Using the list of likely children, for each parent halo figure
        # out the fraction of its particles in the children.
        
        for i,halo in enumerate(parents):
            parents[i]['neighbors_frac'] = {}
            #for neighbor in halo['neighbors']:
            print i,halo['neighbors']
            for j,child in enumerate(children):
                frac = float(len(halo['IDs'] & child['IDs'])) / len(halo['IDs'])
                if frac > 0.3:
                    print 'parent halo ', i, ' has given ', frac,' to child ', j
        return
        
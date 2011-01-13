"""
A very simple, purely-serial, merger tree script that knows how to parse FOF
catalogs output by Enzo and then compare parent/child relationships.

Author: Matthew J. Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

# First pass at a simplified merger tree
#
# Basic outline:
#
# 1. Halo find inline, obtaining particle catalogs
# 2. Load dataset at time t
# 3. Load dataset at time t+1
# 4. Parse catalogs for t and t+1
# 5. Place halos for t+1 in kD-tree
# 6. For every halo in t, execute ball-query with some linking length
# 7. For every halo in ball-query result, execute numpy's intersect1d on
#    particle IDs
# 8. Parentage is described by a fraction of particles that pass from one to
#    the other; we have both descendent fractions and ancestory fractions. 

import numpy as na
import h5py
import time
import pdb
import cPickle

from yt.funcs import *
from yt.utilities.pykdtree import KDTree

# We don't currently use this, but we may again find a use for it in the
# future.
class MaxLengthDict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.order = [None] * 50

    def __setitem__(self, key, val):
        if key not in self.order:
            to_remove = self.order.pop(0)
            self.pop(to_remove, None)
        self.order.append(key)
        dict.__setitem__(self, key, val)

    def __getitem__(self, key):
        if key in self.order:
            self.order.pop(self.order.index(key))
            self.order.append(key)
        return dict.__getitem__(self, key)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self.order.pop(self.order.index(key))
        self.order.insert(0, None)

class HaloCatalog(object):
    cache = None
    def __init__(self, output_id, cache = True):
        r"""A catalog of halos, parsed from EnzoFOF outputs.

        This class will read in catalogs output by the Enzo FOF halo finder and
        make available their positions, radii, etc.  Enzo FOF was provided
        starting with 2.0, and can be run either inline (with the correct
        options) or as a postprocessing step using the `-F` command line
        option.  This class is mostly useful when calculating a merger tree,
        and when the particle IDs for members of a given halo are output as
        well.

        Parameters
        ----------
        output_id : int
            This is the integer output id of the halo catalog to parse and
            load.
        cache : bool
            Should we store, in between accesses, the particle IDs?  If set to
            true, the correct particle files must exist.
        """
        self.output_id = output_id
        self.particle_file = h5py.File("FOF/particles_%04i.h5" % output_id, "r")
        self.parse_halo_catalog()
        if cache: self.cache = dict()#MaxLengthDict()

    def parse_halo_catalog(self):
        hp = []
        for line in open("FOF/groups_%04i.dat" % self.output_id):
            if line.strip() == "": continue # empty
            if line[0] == "#": continue # comment
            if line[0] == "d": continue # datavar
            x,y,z = [float(f) for f in line.split(None, 3)[:-1]]
            hp.append([x,y,z])
        self.halo_positions = na.array(hp)
        self.halo_kdtree = KDTree(self.halo_positions)
        return hp

    def read_particle_ids(self, halo_id):
        if self.cache is not None:
            if halo_id not in self.cache:
                self.cache[halo_id] = self.particle_file["/Halo%08i/Particle ID" % halo_id][:]
            ids = self.cache[halo_id]
        else:
            ids = self.particle_file["/Halo%08i/Particle ID" % halo_id][:]
        return HaloParticleList(halo_id, self.halo_positions[halo_id,:], ids)

    def calculate_parentage_fractions(self, other_catalog, radius = 0.10):
        parentage_fractions = {}
        mylog.debug("Ball-tree query with radius %0.3e", radius)
        all_nearest = self.halo_kdtree.query_ball_tree(
            other_catalog.halo_kdtree, radius)
        pbar = get_pbar("Halo Mergers", HC1.halo_positions.shape[0])
        for hid1, nearest in enumerate(all_nearest):
            pbar.update(hid1)
            parentage_fractions[hid1] = {}
            HPL1 = self.read_particle_ids(hid1)
            for hid2 in nearest:
                HPL2 = other_catalog.read_particle_ids(hid2)
                p1, p2 = HPL1.find_relative_parentage(HPL2)
                parentage_fractions[hid1][hid2] = (p1, p2)
        pbar.finish()
        return parentage_fractions

class HaloParticleList(object):
    def __init__(self, halo_id, position, particle_ids):
        self.halo_id = halo_id
        self.position = na.array(position)
        self.particle_ids = particle_ids

    def find_nearest(self, other_tree, radius = 0.10):
        return other_tree.query_ball_point(self.position, radius)

    def find_relative_parentage(self, child):
        # Return two values: percent this halo gave to the other, and percent
        # of the other that comes from this halo
        overlap = na.intersect1d(self.particle_ids, child.particle_ids).size
        of_child_from_me = float(overlap)/child.particle_ids.size
        of_mine_from_me = float(overlap)/self.particle_ids.size
        return of_child_from_me, of_mine_from_me

def find_halo_relationships(output1_id, output2_id, output_basename = None,
                            radius = 0.10):
    r"""Calculate the parentage and child relationships between two EnzoFOF
    halo catalogs.

    This function performs a very simple merger tree calculation between two
    sets of halos.  For every halo in the second halo catalog, it looks to the
    first halo catalog to find the parents by looking at particle IDs.  The
    particle IDs from the child halos are identified in potential parents, and
    then both percent-of-parent and percent-to-child values are recorded.

    Note that this works only with catalogs constructed by Enzo's FOF halo
    finder, not with catalogs constructed by yt.

    Parameters
    ----------
    output1_id : int
        This is the integer output id of the (first) halo catalog to parse and
        load.
    output2_id : int
        This is the integer output id of the (second) halo catalog to parse and
        load.
    output_basename : string
        If provided, both .cpkl and .txt files containing the parentage
        relationships will be output.
    radius : float, default to 0.10
        In absolute units, the radius to examine when guessing possible
        parent/child relationships.  If this value is too small, you will miss
        possible relationships.

    Returns
    -------
    pfrac : dict
        This is a dict of dicts.  The first key is the parent halo id, the
        second is the child halo id.  The values are the percent contributed
        from parent to child and the percent of a child that came from the
        parent.
    """
    mylog.info("Parsing Halo Catalog %04i", output1_id)
    HC1 = HaloCatalog(output1_id, False)
    mylog.info("Parsing Halo Catalog %04i", output2_id)
    HC2 = HaloCatalog(output2_id, True)
    mylog.info("Calculating fractions")
    pfrac = HC1.calculate_parentage_fractions(HC2)

    if output_basename is not None:
        f = open("%s.txt" % (output_basename), "w")
        for hid1 in sorted(pfrac):
            for hid2 in sorted(pfrac[hid1]):
                p1, p2 = pfrac[hid1][hid2]
                if p1 == 0.0: continue
                f.write( "Halo %s (%s) contributed %0.3e of its particles to %s (%s), which makes up %0.3e of that halo\n" % (
                            hid1, output1_id, p2, hid2, output2_id, p1))
        f.close()

        cPickle.dump(pfrac, open("%s.cpkl" % (output_basename), "wb"))

    return pfrac

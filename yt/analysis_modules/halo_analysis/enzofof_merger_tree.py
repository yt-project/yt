"""
A very simple, purely-serial, merger tree script that knows how to parse FOF
catalogs, either output by Enzo or output by yt's FOF halo finder, and then 
compare parent/child relationships.



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2020, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
# plot_halo_evolution() gives a good full example of how to use the framework

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
#    the other; we have both descendent fractions and ancestor fractions. 


import numpy as np
from yt.utilities.on_demand_imports import _h5py as h5py
import glob
import os

from yt.extern.six.moves import cPickle
from yt.extern.pykdtree import KDTree
from yt.funcs import mylog, get_pbar

import yt.extern.pydot as pydot

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
    external_FOF : bool, optional
        Are we building a tree from outputs generated by an
        external FOF program, or an FOF internal to yt?
    FOF_directory : str, optional
        Directory where FOF files are located
    """
    cache = None
    def __init__(self, output_id, cache = True, external_FOF=True, FOF_directory="FOF"):
        self.output_id = output_id
        self.external_FOF = external_FOF
        self.redshift = 0.0
        self.FOF_directory = FOF_directory
        self.particle_file = h5py.File("%s/particles_%05i.h5" % \
                                      (FOF_directory, output_id), "r")
        if self.external_FOF:
            self.parse_halo_catalog_external()
        else:
            self.parse_halo_catalog_internal()
        if cache: self.cache = dict()#MaxLengthDict()

    def __del__(self):
        self.particle_file.close()

    def parse_halo_catalog_external(self):
        hp = []
        for line in open("%s/groups_%05i.dat" % \
                        (self.FOF_directory, self.output_id)):
            if line.strip() == "": continue # empty
            if line.startswith("# Red"):
                self.redshift = float(line.split("=")[1])
            if line[0] == "#": continue # comment
            if line[0] == "d": continue # datavar
            x,y,z = [float(f) for f in line.split(None, 3)[:-1]]
            hp.append([x,y,z])
        if hp != []:
            self.halo_positions = np.array(hp)
            self.halo_kdtree = KDTree(self.halo_positions)
        else:
            self.halo_positions = None
            self.halo_kdtree = None
        return hp

    def parse_halo_catalog_internal(self):
        """
        This parser works on the files output directly out of yt's internal
        halo_finder.  The parse_halo_catalog_external works with an 
        external version of FOF.

        Examples
        --------
        >>> ds = load("DD0000/DD0000")
        >>> halo_list = FOFHaloFinder(ds)
        >>> halo_list.write_out("FOF/groups_00000.txt")
        >>> halos_COM = parse_halo_catalog_internal()
        """
        hp = []
        for line in open("%s/groups_%05i.txt" % \
                        (self.FOF_directory, self.output_id)):
            if line.startswith("# RED"):
                self.redshift = float(line.split("=")[1])
                continue
            if line.strip() == "": continue # empty
            if line[0] == "#": continue # comment
            x,y,z = [float(f) for f in line.split()[7:10]] # COM x,y,z
            hp.append([x,y,z])
        if hp != []:
            self.halo_positions = np.array(hp)
            self.halo_kdtree = KDTree(self.halo_positions)
        else:
            self.halo_positions = None
            self.halo_kdtree = None
        return hp

    def read_particle_ids(self, halo_id):        
        if self.cache is not None:
            if halo_id not in self.cache:
                if self.external_FOF:
                    self.cache[halo_id] = \
                    self.particle_file["/Halo%08i/Particle ID" % halo_id][:]
                else:
                    self.cache[halo_id] = \
                    self.particle_file["/Halo%08i/particle_index" % halo_id][:]
            ids = self.cache[halo_id]
        else:
            if self.external_FOF:
                ids = self.particle_file["/Halo%08i/Particle ID" % halo_id][:]
            else:
                ids = self.particle_file["/Halo%08i/particle_index" % halo_id][:]
        return HaloParticleList(halo_id, self.halo_positions[halo_id,:], ids)

    def calculate_parentage_fractions(self, other_catalog, radius = 0.10):
        parentage_fractions = {}
        if self.halo_positions is None or other_catalog.halo_positions is None:
            return parentage_fractions
        mylog.debug("Ball-tree query with radius %0.3e", radius)
        all_nearest = self.halo_kdtree.query_ball_tree(
            other_catalog.halo_kdtree, radius)
        pbar = get_pbar("Halo Mergers", self.halo_positions.shape[0])
        for hid1, nearest in enumerate(all_nearest):
            pbar.update(hid1)
            parentage_fractions[hid1] = {}
            HPL1 = self.read_particle_ids(hid1)
            for hid2 in sorted(nearest):
                HPL2 = other_catalog.read_particle_ids(hid2)
                p1, p2 = HPL1.find_relative_parentage(HPL2)
                parentage_fractions[hid1][hid2] = (p1, p2, HPL2.number_of_particles)
            parentage_fractions[hid1]["NumberOfParticles"] = HPL1.number_of_particles
        pbar.finish()
        return parentage_fractions

class HaloParticleList(object):
    def __init__(self, halo_id, position, particle_ids):
        self.halo_id = halo_id
        self.position = np.array(position)
        self.particle_ids = particle_ids
        self.number_of_particles = particle_ids.size

    def find_nearest(self, other_tree, radius = 0.10):
        return other_tree.query_ball_point(self.position, radius)

    def find_relative_parentage(self, child):
        # Return two values: percent this halo gave to the other, and percent
        # of the other that comes from this halo
        overlap = np.intersect1d(self.particle_ids, child.particle_ids).size
        of_child_from_me = float(overlap)/child.particle_ids.size
        of_mine_from_me = float(overlap)/self.particle_ids.size
        return of_child_from_me, of_mine_from_me

class EnzoFOFMergerBranch(object):
    def __init__(self, tree, output_num, halo_id, max_children,
                 min_relation=0.25):
        self.output_num = output_num
        self.halo_id = halo_id
        self.npart = tree.relationships[output_num][halo_id]["NumberOfParticles"]
        self.children = []
        self.progenitor = -1
        max_relationship = 0.0
        halo_count = 0
        keys = list(tree.relationships[output_num][halo_id].keys())
        keys.remove('NumberOfParticles')
        for k in sorted(keys):
            v = tree.relationships[output_num][halo_id][k]
            if v[1] > min_relation and halo_count < max_children:
                halo_count += 1
                self.children.append((k,v[1],v[2]))
                if v[1] > max_relationship:
                    self.progenitor = k
                    max_relationship = v[1]

class EnzoFOFMergerTree(object):
    r"""Calculates the parentage relationships for halos for a series of
    outputs, using the framework provided in enzofof_merger_tree.

    Parameters
    ----------
    zrange : tuple
        This is the redshift range (min, max) to calculate the
        merger tree. E.g. (0, 2) for z=2 to z=0
    cycle_range : tuple, optional
        This is the cycle number range (min, max) to calculate the
        merger tree.  If both zrange and cycle_number given,
        ignore zrange.
    output : bool, optional
        If provided, both .cpkl and .txt files containing the parentage
        relationships will be output.
    load_saved : bool, optional
        Flag to load previously saved parental relationships
    save_filename : str, optional
        Filename to save parental relationships
    external_FOF : bool, optional
        Are we building a tree from outputs generated by an
        external FOF program, or an FOF internal to yt?
    FOF_directory : str, optional
        Directory where FOF files are located, note that the files
        must be named according to the syntax: groups_DDDDD.txt for
        internal yt outputs, and groups_DDDDD.dat for external FOF outputs.
        where DDDDD are digits representing the equivalent cycle number.
        e.g. groups_00000.txt
    
    Examples
    --------
    >>> mt = EnzoFOFMergerTree()    # by default it grabs every DD in FOF dir
    >>> mt.build_tree(0)  # Create tree for halo 0
    >>> mt.print_tree()
    >>> mt.write_dot()

    See Also
    --------
    plot_halo_evolution()
    """    
    def __init__(self, zrange=None, cycle_range=None, output=False,
                 load_saved=False, save_filename="merger_tree.cpkl",
                 external_FOF=True, FOF_directory="FOF"):

        self.relationships = {}
        self.redshifts = {}
        self.external_FOF = external_FOF
        self.FOF_directory = FOF_directory
        if load_saved:
            self.load_tree("%s/%s" % (self.FOF_directory, save_filename))
            # make merger tree work within specified cycle/z limits
            # on preloaded halos
            if zrange is not None:
                self.select_redshifts(zrange)
            if cycle_range is not None:
                self.select_cycles(cycle_range)
        else:
            self.find_outputs(zrange, cycle_range, output)
            self.run_merger_tree(output)
            self.save_tree("%s/%s" % (self.FOF_directory, save_filename))
        
    def select_cycles(self, cycle_range):
        """
        Takes an existing tree and pares it to only include a subset of
        cycles.  Useful in paring a loaded tree. 
        """
        # N.B. Does not delete info from self.relationships to save space
        # just removes it from redshift dict for indexing
        for cycle in self.redshifts.keys():
            if cycle <= cycle_range[0] and cycle >= cycle_range[1]:
                del self.redshifts[cycle]

    def select_redshifts(self, zrange):
        """
        Takes an existing tree and pares it to only include a subset of
        redshifts.  Useful in paring a loaded tree. 
        """
        # N.B. Does not delete info from self.relationships to save space
        # just removes it from redshift dict for indexing
        for redshift in self.redshifts.values():
            if redshift <= zrange[0] and redshift >= zrange[1]:
                # some reverse lookup magic--assumes unique cycle/z pairs
                cycle = [key for key,value in self.redshifts.items() \
                         if value == redshift][0]
                del self.redshifts[cycle]

    def save_tree(self, filename):
        cPickle.dump((self.redshifts, self.relationships),
                     open(filename, "wb"))

    def load_tree(self, filename):
        self.redshifts, self.relationships = \
                        cPickle.load(open(filename, "rb"))

    def clear_data(self):
        r"""Deletes previous merger tree, but keeps parentage
        relationships.
        """
        del self.levels

    def find_outputs(self, zrange, cycle_range, output):
        self.numbers = []
        if self.external_FOF:
            filenames = "%s/groups_*.dat" % (self.FOF_directory)
            files = glob.glob(filenames)
        else:
            filenames = "%s/groups_*.txt" % (self.FOF_directory)
            files = glob.glob(filenames)
        # If using redshift range, load redshifts only
        for f in files:
            num = int(f[-9:-4])
            if zrange is not None:
                HC = HaloCatalog(num, external_FOF=self.external_FOF, \
                                FOF_directory=self.FOF_directory)
                # Allow for some epsilon
                diff1 = (HC.redshift - zrange[0]) / zrange[0]
                diff2 = (HC.redshift - zrange[1]) / zrange[1]
                if diff1 >= -1e-3 and diff2 <= 1e-3:
                    self.numbers.append(num)
                del HC
            elif cycle_range is not None:
                if num >= cycle_range[0] and num <= cycle_range[1]:
                    self.numbers.append(num)
            else:
                self.numbers.append(num)
        self.numbers.sort()

    def run_merger_tree(self, output):
        # Run merger tree for all outputs, starting with the last output
        for i in range(len(self.numbers)-1, 0, -1):
            if output:
                output = "%s/tree-%5.5d-%5.5d" % \
                         (self.FOF_directory, self.numbers[i], self.numbers[i-1])
            else:
                output = None
            z0, z1, fr = find_halo_relationships(self.numbers[i], \
                                                 self.numbers[i-1], \
                                                 output_basename=output, \
                                                 external_FOF=self.external_FOF,
                                                 FOF_directory=self.FOF_directory)
            self.relationships[self.numbers[i]] = fr
            self.redshifts[self.numbers[i]] = z0
        # Fill in last redshift
        self.redshifts[self.numbers[0]] = z1

    def build_tree(self, halonum, min_particles=0, max_children=1e20):
        r"""Builds a merger tree, starting at the last output.

        Parameters
        ----------
        halonum : int
            Halo number in the last output to analyze.
        min_particles : int, optional
            Minimum number of particles of halos in tree.
        max_children : int, optional
            Maximum number of child halos each leaf can have.
        """
        self.halonum = halonum
        self.max_children = max_children
        self.output_numbers = sorted(self.relationships, reverse=True)
        self.levels = {}
        trunk = self.output_numbers[0]
        self.levels[trunk] = [EnzoFOFMergerBranch(self, trunk, halonum,
                                                  max_children)]
        self.generate_tree(min_particles, max_children)

    def filter_small_halos(self, lvl, min_particles):
        # Filter out children with less than min_particles
        for h in self.levels[lvl]:
            fil = []
            for c in h.children:
                if c[2] > min_particles:  # c[2] = npart
                    fil.append(c)
            h.children = fil

    def generate_tree(self, min_particles, max_children):
        self.filter_small_halos(self.output_numbers[0], min_particles)
        for i in range(1,len(self.output_numbers)):
            prev = self.output_numbers[i-1]
            this = self.output_numbers[i]
            self.levels[this] = []
            this_halos = []  # To check for duplicates
            for h in self.levels[prev]:
                for c in h.children:
                    if c[0] in this_halos: continue
                    if self.relationships[this] == {}: continue
                    branch = EnzoFOFMergerBranch(self, this, c[0],
                                                 max_children)
                    self.levels[this].append(branch)
                    this_halos.append(c[0])
            self.filter_small_halos(this, min_particles)

    def get_massive_progenitors(self, halonum, min_relation=0.25):
        r"""Returns a list of the most massive progenitor halos.

        This routine walks down the tree, following the most massive
        progenitor on each node.

        Parameters
        ----------
        halonum : int
            Halo number at the last output to trace.

        Returns
        -------
        output : dict
            Dictionary of redshifts, cycle numbers, and halo numbers
            of the most massive progenitor.  keys = {redshift, cycle,
            halonum}
        """
        output = {"redshift": [], "cycle": [], "halonum": []}
        # First (lowest redshift) node in tree
        halo0 = halonum
        for cycle in sorted(self.numbers, reverse=True):
            if cycle not in self.relationships: break
            if halo0 not in self.relationships[cycle]: break
            node = self.relationships[cycle][halo0]
            output["redshift"].append(self.redshifts[cycle])
            output["cycle"].append(cycle)
            output["halonum"].append(halo0)
            # Find progenitor
            max_rel = 0.0
            for k,v in node.items():
                if not str(k).isdigit(): continue
                if v[1] > max_rel and v[1] > min_relation:
                    halo0 = k
                    max_rel = v[1]
        return output

    def print_tree(self):
        r"""Prints the merger tree to stdout.
        """
        for lvl in sorted(self.levels, reverse=True):
            if lvl not in self.redshifts: continue
            print("========== Cycle %5.5d (z=%f) ==========" % \
                  (lvl, self.redshifts[lvl]))
            for br in self.levels[lvl]:
                print("Parent halo = %d" % br.halo_id)
                print("--> Most massive progenitor == Halo %d" % \
                      (br.progenitor))
                for i,c in enumerate(br.children):
                    if i > self.max_children: break
                    print("-->    Halo %8.8d :: fraction = %g" % (c[0], c[1]))

    def save_halo_evolution(self, filename):
        """
        Saves as an HDF5 file the relevant details about a halo
        over the course of its evolution following the most massive
        progenitor to have given it the bulk of its particles.
        It stores info from the FOF_groups file: location, mass, id, etc.
        """
        f = h5py.File("%s/%s" % (self.FOF_directory, filename), 'a')
        cycle_fin = sorted(list(self.redshifts.keys()))[-1]
        halo_id = self.levels[cycle_fin][0].halo_id
        halo = "halo%05d" % halo_id
        if halo in f:
            del f["halo%05d" % halo_id]
        g = f.create_group("halo%05d" % halo_id)
        size = len(self.redshifts)
        cycle = np.zeros(size)
        redshift = np.zeros(size)
        halo_id = np.zeros(size)
        fraction = np.zeros(size)
        mass = np.zeros(size)
        densest_point = np.zeros((3,size))
        COM = np.zeros((6,size))
        fraction[0] = 1.

        for i, lvl in enumerate(sorted(self.levels, reverse=True)):
            if len(self.levels[lvl]) == 0:  # lineage for this halo ends
                cycle = cycle[:i]           # so truncate arrays, and break
                redshift = redshift[:i]     # Not big enough.
                halo_id = halo_id[:i]
                fraction = fraction[:i]
                mass = mass[:i]
                densest_point = densest_point[:,:i]
                COM = COM[:,:i]
                break   
            if lvl not in self.redshifts: continue
            mylog.info("========== Cycle %5.5d (z=%f) ==========" % \
                  (lvl, self.redshifts[lvl]))
            cycle[i] = lvl 
            redshift[i] = self.redshifts[lvl]

            br = self.levels[lvl][0]
            mylog.info("Parent halo = %d" % br.halo_id)
            mylog.info("-> Most massive progenitor == Halo %d" % (br.progenitor))
            halo_id[i] = br.halo_id

            if len(br.children) == 0:     # lineage for this halo ends 
                cycle = cycle[:i+1]       # (no children)
                redshift = redshift[:i+1] # so truncate arrays, and break
                halo_id = halo_id[:i+1]
                fraction = fraction[:i+1]
                mass = mass[:i+1]
                densest_point = densest_point[:,:i+1]
                COM = COM[:,:i+1]
                break   

            if i < size-1:
                fraction[i+1] = br.children[0][1]  

            # open up FOF file to parse for details
            filename = "%s/groups_%05d.txt" % (self.FOF_directory, lvl)
            mass[i], densest_point[:,i], COM[:,i] = \
                    grab_FOF_halo_info_internal(filename, br.halo_id)

        # save the arrays in the hdf5 file
        g.create_dataset("cycle", data=cycle)
        g.create_dataset("redshift", data=redshift)
        g.create_dataset("halo_id", data=halo_id)
        g.create_dataset("fraction", data=fraction)
        g.create_dataset("mass", data=mass)
        g.create_dataset("densest_point", data=densest_point)
        g.create_dataset("COM", data=COM)
        f.close()

    def write_dot(self, filename=None):
        r"""Writes merger tree to a GraphViz or image file.

        Parameters
        ----------
        filename : str, optional
            Filename to write the GraphViz file.  Default will be
            tree_halo%05i.gv, which is a text file in the GraphViz format.
            If filename is an image (e.g. "MergerTree.png") the output will
            be in the appropriate image format made by calling GraphViz
            automatically. See GraphViz (e.g. "dot -v")
            for a list of available output formats.
        """
        if filename is None:
            filename = "%s/tree_halo%5.5d.gv" % \
                        (self.FOF_directory, self.halonum)
        # Create the pydot graph object.
        self.graph = pydot.Dot('galaxy', graph_type='digraph')
        self.halo_shape = "rect"
        self.z_shape = "plaintext"
        # Subgraphs to align levels
        self.subgs = {}
        for num in self.numbers:
            self.subgs[num] = pydot.Subgraph('', rank = 'same')
            self.graph.add_subgraph(self.subgs[num])
        sorted_lvl = sorted(self.levels, reverse=True)
        for ii,lvl in enumerate(sorted_lvl):
            # Since we get the cycle number from the key, it won't
            # exist for the last level, i.e. children of last level.
            # Get it from self.numbers.
            if ii < len(sorted_lvl)-1:
                next_lvl = sorted_lvl[ii+1]
            else:
                next_lvl = self.numbers[0]
            for br in self.levels[lvl]:
                for c in br.children:
                    color = "red" if c[0] == br.progenitor else "black"
                    self.graph.add_edge(pydot.Edge("C%d_H%d" %(lvl, br.halo_id),
                        "C%d_H%d" % (next_lvl, c[0]), color=color))
                    #line = "    C%d_H%d -> C%d_H%d [color=%s];\n" % \
                    #      (lvl, br.halo_id, next_lvl, c[0], color)
                    
                    #fp.write(line)
        for ii,lvl in enumerate(sorted_lvl):
            npart_max = 0
            for br in self.levels[lvl]:
                if br.npart > npart_max: npart_max = br.npart
            for br in self.levels[lvl]:
                halo_str = "C%d_H%d" % (lvl, br.halo_id)
                style = "filled" if br.npart == npart_max else "solid"
                self.graph.add_node(pydot.Node(halo_str,
                label = "Halo %d\\n%d particles" % (br.halo_id, br.npart),
                style = style, shape = self.halo_shape))
                # Add this node to the correct level subgraph.
                self.subgs[lvl].add_node(pydot.Node(halo_str))
        for lvl in self.numbers:
            # Don't add the z if there are no halos already in the subgraph.
            if len(self.subgs[lvl].get_node_list()) == 0: continue
            self.subgs[lvl].add_node(pydot.Node("%1.5e" % self.redshifts[lvl],
                shape = self.z_shape, label = "z=%0.3f" % self.redshifts[lvl]))
        # Based on the suffix of the file name, write out the result to a file.
        suffix = filename.split(".")[-1]
        if suffix == "gv": suffix = "raw"
        mylog.info("Writing %s format %s to disk." % (suffix, filename))
        self.graph.write("%s" % filename, format=suffix)

def find_halo_relationships(output1_id, output2_id, output_basename = None,
                            radius = 0.10, external_FOF=True, 
                            FOF_directory='FOF'):
    r"""Calculate the parentage and child relationships between two EnzoFOF
    halo catalogs.

    This function performs a very simple merger tree calculation between two
    sets of halos.  For every halo in the second halo catalog, it looks to the
    first halo catalog to find the parents by looking at particle IDs.  The
    particle IDs from the child halos are identified in potential parents, and
    then both percent-of-parent and percent-to-child values are recorded.

    Note that this works with catalogs constructed by Enzo's FOF halo
    when used in external_FOF=True mode, whereas it will work with 
    catalogs constructed by yt using external_FOF=False mode.

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
    FOF_directory : str, optional
        Directory where FOF files are located

    Returns
    -------
    pfrac : dict
        This is a dict of dicts.  The first key is the parent halo id, the
        second is the child halo id.  The values are the percent contributed
        from parent to child and the percent of a child that came from the
        parent.
    """
    mylog.info("Parsing Halo Catalog %04i", output1_id)
    HC1 = HaloCatalog(output1_id, False, external_FOF=external_FOF, \
                      FOF_directory=FOF_directory)
    mylog.info("Parsing Halo Catalog %04i", output2_id)
    HC2 = HaloCatalog(output2_id, True, external_FOF=external_FOF, \
                      FOF_directory=FOF_directory)
    mylog.info("Calculating fractions")
    pfrac = HC1.calculate_parentage_fractions(HC2)

    if output_basename is not None and pfrac != {}:
        f = open("%s.txt" % (output_basename), "w")
        for hid1 in sorted(pfrac):
            for hid2 in sorted(pfrac[hid1]):
                if not str(hid2).isdigit(): continue
                p1, p2, npart = pfrac[hid1][hid2]
                if p1 == 0.0: continue
                f.write( "Halo %s (%s) contributed %0.3e of its particles to %s (%s), which makes up %0.3e of that halo\n" % (
                            hid1, output1_id, p2, hid2, output2_id, p1))
        f.close()

        cPickle.dump(pfrac, open("%s.cpkl" % (output_basename), "wb"))

    return HC1.redshift, HC2.redshift, pfrac

def grab_FOF_halo_info_internal(filename, halo_id):
    """
    Finds a specific halo's information in the FOF group output information
    and pass relevant parameters to caller.
    """
    # open up FOF file to parse for details
    groups_file = open(filename, 'r')
    for line in groups_file:
        if line.startswith("#"): continue
        if int(line.split()[0]) == halo_id:
            ar = np.array(line.split()).astype('float64')
            return ar[1], ar[4:7], ar[7:13]  # mass, xyz_dens, xyzvxvyvz_COM

def plot_halo_evolution(filename, halo_id, x_quantity='cycle', y_quantity='mass',
                        x_log=False, y_log=True, FOF_directory='FOF'):
    """
    Once you have generated a file using the 
    EnzoFOFMergerTree.save_halo_evolution function, this is a simple way of 
    plotting the evolution in the quantities of that halo over its lifetime.

    Parameters
    ----------
    filename : str
        The filename to which you saved the hdf5 data from save_halo_evolution
    halo_id : int
        The halo in 'filename' that you want to follow
    x_quantity : str, optional
        The quantity that you want to plot as the x_coord.
        Valid options are:

           * cycle
           * mass
           * fraction
           * halo_id
           * redshift
           * dense_x
           * dense_y
           * dense_z
           * COM_x
           * COM_y
           * COM_z
           * COM_vx
           * COM_vy
           * COM_vz

    y_quantity : str, optional
        The quantity that you want to plot as the y_coord.
    x_log : bool, optional
        Do you want the x-axis to be in log or linear?
    y_log : bool, optional
        Do you want the y-axis to be in log or linear?
    FOF_directory : str, optional
        Directory where FOF files (and hdf file) are located

    Examples
    --------

    >>> # generates mass history plots for the 20 most massive halos at t_fin.
    >>> ts = DatasetSeries.from_filenames("DD????/DD????")
    >>> # long step--must run FOF on each DD, but saves outputs for later use
    >>> for ds in ts:   
    ...     halo_list = FOFHaloFinder(ds)
    ...     i = int(ds.basename[2:])
    ...     halo_list.write_out("FOF/groups_%05i.txt" % i)
    ...     halo_list.write_particle_lists("FOF/particles_%05i" % i)
    ...
    >>> mt = EnzoFOFMergerTree(external_FOF=False)
    >>> for i in range(20):
    ...     mt.build_tree(i)
    ...     mt.save_halo_evolution('halos.h5')
    ...
    >>> for i in range(20):
    ...     plot_halo_evolution('halos.h5', i)
    """
    import matplotlib.pyplot as plt
    f = h5py.File("%s/%s" % (FOF_directory, filename), 'r')
    basename = os.path.splitext(filename)[0]
    halo = "halo%05d" % halo_id
    basename = basename + "_" + halo
    g = f[halo]
    values = list(g)
    index_dict = {'x' : 0, 'y' : 1, 'z' : 2, 'vx' : 3, 'vy' : 4, 'vz' : 5}
    coords = {}
    fields = {}
    for i, quantity in enumerate((x_quantity, y_quantity)):
        field = quantity
        if quantity.startswith('COM'):
            index = index_dict[quantity.split('_')[-1]]
            quantity = ('COM')
        if quantity.startswith('dense'):
            index = index_dict[quantity.split('_')[-1]]
            quantity = ('densest_point')
        if quantity not in values:
            exit('%s not in list of values in %s for halo %d' % \
                 (quantity, filename, halo_id))
        if not field == quantity:
            coords[i] = g[quantity][index,:]
        else:
            coords[i] = g[quantity]
        if len(coords[i]) == 1: 
            # ("Only 1 value for Halo %d.  Ignoring." % halo_id)
            return
        fields[i] = field

    ax = plt.axes()
    ax.plot(coords[0], coords[1])
    ax.set_title(basename)
    ax.set_xlabel(fields[0])
    ax.set_ylabel(fields[1])
    if x_log:
        ax.set_xscale("log")
    if y_log:
        ax.set_yscale("log")
    ofn = "%s/%s_%s_%s.png" % (FOF_directory, basename, fields[0], fields[1])
    plt.savefig(ofn)
    plt.clf()

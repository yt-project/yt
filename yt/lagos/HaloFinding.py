"""
HOP-output data handling

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Stephen Skory <stephenskory@yahoo.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *
from yt.lagos.hop.EnzoHop import RunHOP
from yt.lagos.chainHOP.chainHOP import *
try:
    from yt.lagos.fof.EnzoFOF import RunFOF
except ImportError:
    pass

from kd import *
import math

class Halo(object):
    """
    A data source that returns particle information about the members of a
    HOP-identified halo.
    """
    __metaclass__ = ParallelDummy # This will proxy up our methods
    _distributed = False
    _processing = False
    _owner = 0
    indices = None
    dont_wrap = ["get_sphere", "write_particle_list"]
    extra_wrap = ["__getitem__"]

    def __init__(self, halo_list, id, indices = None):
        self.halo_list = halo_list
        self.id = id
        self.data = halo_list._data_source
        if indices is not None: self.indices = halo_list._base_indices[indices]
        # We assume that if indices = None, the instantiator has OTHER plans
        # for us -- i.e., setting it somehow else

    def center_of_mass(self):
        """
        Calculate and return the center of mass.
        """
        c_vec = self.maximum_density_location() - na.array([0.5,0.5,0.5])
        pm = self["ParticleMassMsun"]
        cx = (self["particle_position_x"] - c_vec[0])
        cy = (self["particle_position_y"] - c_vec[1])
        cz = (self["particle_position_z"] - c_vec[2])
        com = na.array([v-na.floor(v) for v in [cx,cy,cz]])
        return (com*pm).sum(axis=1)/pm.sum() + c_vec

    def maximum_density(self):
        """
        Return the HOP-identified maximum density.
        """
        return self.halo_list._max_dens[self.id][0]

    def maximum_density_location(self):
        """
        Return the location HOP identified as maximally dense.
        """
        return na.array([
                self.halo_list._max_dens[self.id][1],
                self.halo_list._max_dens[self.id][2],
                self.halo_list._max_dens[self.id][3]])

    def total_mass(self):
        """
        Returns the total mass in solar masses of the halo.
        """
        return self["ParticleMassMsun"].sum()

    def bulk_velocity(self):
        """
        Returns the mass-weighted average velocity.
        """
        pm = self["ParticleMassMsun"]
        vx = (self["particle_velocity_x"] * pm).sum()
        vy = (self["particle_velocity_y"] * pm).sum()
        vz = (self["particle_velocity_z"] * pm).sum()
        return na.array([vx,vy,vz])/pm.sum()

    def maximum_radius(self, center_of_mass=True):
        """
        Returns the maximum radius in the halo for all particles,
        either from the point of maximum density or from the (default)
        *center_of_mass*.
        """
        if center_of_mass: center = self.center_of_mass()
        else: center = self.maximum_density_location()
        rx = na.abs(self["particle_position_x"]-center[0])
        ry = na.abs(self["particle_position_y"]-center[1])
        rz = na.abs(self["particle_position_z"]-center[2])
        r = na.sqrt(na.minimum(rx, 1.0-rx)**2.0
                +   na.minimum(ry, 1.0-ry)**2.0
                +   na.minimum(rz, 1.0-rz)**2.0)
        return r.max()

    def __getitem__(self, key):
        return self.data[key][self.indices]

    def get_sphere(self, center_of_mass=True):
        """
        Returns an EnzoSphere centered on either the point of maximum density
        or the *center_of_mass*, with the maximum radius of the halo.
        """
        if center_of_mass: center = self.center_of_mass()
        else: center = self.maximum_density_location()
        radius = self.maximum_radius()
        # A bit of a long-reach here...
        sphere = self.halo_list._data_source.hierarchy.sphere(
                        center, radius=radius)
        return sphere

    def get_size(self):
        return self.indices.size

    def write_particle_list(self, handle):
        self._processing = True
        gn = "Halo%08i" % (self.id)
        handle.create_group("/%s" % gn)
        for field in ["particle_position_%s" % ax for ax in 'xyz'] \
                   + ["particle_velocity_%s" % ax for ax in 'xyz'] \
                   + ["particle_index"]:
            handle.create_dataset("/%s/%s" % (gn, field), data=self[field])
        n = handle["/%s" % gn]
        # set attributes on n
        self._processing = False

class HOPHalo(Halo):
    pass

class chainHOPHalo(Halo):
    pass

class FOFHalo(Halo):

    def center_of_mass(self):
        """
        Calculate and return the center of mass.
        """
        pm = self["ParticleMassMsun"]
        cx = self["particle_position_x"]
        cy = self["particle_position_y"]
        cz = self["particle_position_z"]
        c_vec = na.array([cx[0],cy[0],cz[0]]) - na.array([0.5,0.5,0.5])
        cx = cx - c_vec[0]
        cy = cy - c_vec[1]
        cz = cz - c_vec[2]
        com = na.array([v-na.floor(v) for v in [cx,cy,cz]])
        com = (pm * com).sum(axis=1)/pm.sum() + c_vec
        return com

    def maximum_density(self):
        return -1

    def maximum_density_location(self):
        return self.center_of_mass()

class HaloList(object):

    _fields = ["particle_position_%s" % ax for ax in 'xyz']

    def __init__(self, data_source, dm_only = True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self._data_source = data_source
        self.dm_only = dm_only
        self._groups = []
        self._max_dens = {}
        self.__obtain_particles()
        self._run_finder()
        mylog.info("Parsing outputs")
        self._parse_output()
        mylog.debug("Finished. (%s)", len(self))

    def __obtain_particles(self):
        if self.dm_only: ii = self.__get_dm_indices()
        else: ii = slice(None)
        self.particle_fields = {}
        for field in self._fields:
            tot_part = self._data_source[field].size
            self.particle_fields[field] = self._data_source[field][ii]
        self._base_indices = na.arange(tot_part)[ii]

    def __get_dm_indices(self):
        if 'creation_time' in self._data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on creation time")
            return (self._data_source["creation_time"] < 0)
        elif 'particle_type' in self._data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on particle type")
            return (self._data_source["particle_type"] == 1)
        else:
            mylog.warning("No particle_type, no creation_time, so not distinguishing.")
            return slice(None)

    def _parse_output(self):
        unique_ids = na.unique(self.tags)
        counts = na.bincount(self.tags+1)
        sort_indices = na.argsort(self.tags)
        grab_indices = na.indices(self.tags.shape).ravel()[sort_indices]
        dens = self.densities[sort_indices]
        cp = 0
        for i in unique_ids:
            cp_c = cp + counts[i+1]
            if i == -1:
                cp += counts[i+1]
                continue
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(self._halo_class(self, i, group_indices))
            md_i = na.argmax(dens[cp:cp_c])
            px, py, pz = [self.particle_fields['particle_position_%s'%ax][group_indices]
                                            for ax in 'xyz']
            self._max_dens[i] = (dens[cp:cp_c][md_i], px[md_i], py[md_i], pz[md_i])
            cp += counts[i+1]

    def __len__(self):
        return len(self._groups)
 
    def __iter__(self):
        for i in self._groups: yield i

    def __getitem__(self, key):
        return self._groups[key]

    def nearest_neighbors_3D(self, haloID, num_neighbors=7, search_radius=.2):
        """
        for halo *haloID*, find up to *num_neighbors* nearest neighbors in 3D
        using the kd tree. Search over *search_radius* in code units.
        Returns a list of the neighbors distances and ID with format
        [distance,haloID].
        """
        period = self.pf['DomainRightEdge'] - self.pf['DomainLeftEdge']
        # Initialize the dataset of points from all the haloes
        dataset = []
        for group in self:
            p = Point()
            p.data = group.center_of_mass().tolist()
            p.haloID = group.id
            dataset.append(p)
        mylog.info('Building kd tree...')
        kd = buildKdHyperRectTree(dataset[:],2*num_neighbors)
        # make the neighbors object
        neighbors = Neighbors()
        neighbors.k = num_neighbors
        neighbors.points = []
        neighbors.minDistanceSquared = search_radius * search_radius
        mylog.info('Finding nearest neighbors...')
        getKNN(self[haloID].center_of_mass().tolist(), kd, neighbors,0., period.tolist())
        # convert the data in order to return something less perverse than a
        # Neighbors object, also root the distances
        n_points = []
        for n in neighbors.points:
            n_points.append([math.sqrt(n[0]),n[1].haloID])
        return n_points

    def nearest_neighbors_2D(self, haloID, num_neighbors=7, search_radius=.2,
        proj_dim=0):
        """
        for halo *haloID*, find up to *num_neighbors* nearest neighbors in 2D
        using the kd tree. Search over *search_radius* in code units.
        The halo positions are projected along dimension *proj_dim*.
        Returns a list of the neighbors distances and ID with format
        [distance,haloID].
        """
        # Set up a vector to multiply other vectors by to project along proj_dim
        vec = na.array([1.,1.,1.])
        vec[proj_dim] = 0.
        period = self.pf['DomainRightEdge'] - self.pf['DomainLeftEdge']
        period = period * vec
        # Initialize the dataset of points from all the haloes
        dataset = []
        for group in self:
            p = Point()
            cm = group.center_of_mass() * vec
            p.data = cm.tolist()
            p.haloID = group.id
            dataset.append(p)
        mylog.info('Building kd tree...')
        kd = buildKdHyperRectTree(dataset[:],2*num_neighbors)
        # make the neighbors object
        neighbors = Neighbors()
        neighbors.k = num_neighbors
        neighbors.points = []
        neighbors.minDistanceSquared = search_radius * search_radius
        mylog.info('Finding nearest neighbors...')
        cm = self[haloID].center_of_mass() * vec
        getKNN(cm.tolist(), kd, neighbors,0., period.tolist())
        # convert the data in order to return something less perverse than a
        # Neighbors object, also root the distances
        n_points = []
        for n in neighbors.points:
            n_points.append([math.sqrt(n[0]),n[1].haloID])
        return n_points

    def write_out(self, filename):
        """
        Write out standard HOP information to *filename*.
        """
        if hasattr(filename, 'write'):
            f = filename
        else:
            f = open(filename,"w")
        f.write("# HALOS FOUND WITH %s\n" % (self._name))
        f.write("\t".join(["# Group","Mass","# part","max dens"
                           "x","y","z", "center-of-mass",
                           "x","y","z",
                           "vx","vy","vz","max_r","\n"]))
        for group in self:
            f.write("%10i\t" % group.id)
            f.write("%0.9e\t" % group.total_mass())
            f.write("%10i\t" % group.get_size())
            f.write("%0.9e\t" % group.maximum_density())
            f.write("\t".join(["%0.9e" % v for v in group.maximum_density_location()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.center_of_mass()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.bulk_velocity()]))
            f.write("\t")
            f.write("%0.9e\t" % group.maximum_radius())
            f.write("\n")
            f.flush()
        f.close()

class HOPHaloList(HaloList):

    _name = "HOP"
    _halo_class = HOPHalo
    _fields = ["particle_position_%s" % ax for ax in 'xyz'] + \
              ["ParticleMassMsun"]

    def __init__(self, data_source, threshold=160.0, dm_only=True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self.threshold = threshold
        mylog.info("Initializing HOP")
        HaloList.__init__(self, data_source, dm_only)

    def _run_finder(self):
        self.densities, self.tags = \
            RunHOP(self.particle_fields["particle_position_x"],
                   self.particle_fields["particle_position_y"],
                   self.particle_fields["particle_position_z"],
                   self.particle_fields["ParticleMassMsun"],
                   self.threshold)
        self.particle_fields["densities"] = self.densities
        self.particle_fields["tags"] = self.tags

    def write_out(self, filename="HopAnalysis.out"):
        HaloList.write_out(self, filename)

class FOFHaloList(HaloList):
    _name = "FOF"
    _halo_class = FOFHalo

    def __init__(self, data_source, link=0.2, dm_only=True):
        self.link = link
        mylog.info("Initializing FOF")
        HaloList.__init__(self, data_source, dm_only)

    def _run_finder(self):
        self.tags = \
            RunFOF(self.particle_fields["particle_position_x"],
                   self.particle_fields["particle_position_y"],
                   self.particle_fields["particle_position_z"],
                   self.link)
        self.densities = na.ones(self.tags.size, dtype='float64') * -1
        self.particle_fields["densities"] = self.densities
        self.particle_fields["tags"] = self.tags

    def write_out(self, filename="FOFAnalysis.out"):
        HaloList.write_out(self, filename)

class chainHOPHaloList(HaloList):
    _name = "chainHOP"
    _halo_class = chainHOPHalo
    _fields = ["particle_position_%s" % ax for ax in 'xyz'] + \
              ["ParticleMassMsun"]
    
    def __init__(self, data_source, padding, num_neighbors, bounds, total_mass,
        threshold=160.0, dm_only=True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self.threshold = threshold
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.total_mass = total_mass
        self.period = self._data_source.pf['DomainRightEdge'] - \
            self._data_source.pf['DomainLeftEdge']
        mylog.info("Initializing HOP")
        HaloList.__init__(self, data_source, dm_only)

    def _is_inside(self,point):
        # test to see if this point is inside the 'real' region
        (LE,RE) = self.bounds
        point = na.array(point)
        if (point >= LE).all() and (point < RE).all():
            return True
        else:
            return False

    def _run_finder(self):
        xt = self.particle_fields["particle_position_x"] #- 0.961
        yt = self.particle_fields["particle_position_y"] #- 0.369
        zt = self.particle_fields["particle_position_z"] #- 0.710
        for i,x in enumerate(xt):
            if x < 0:
                xt[i] = 1+x
        for i,y in enumerate(yt):
            if y < 0:
                yt[i] = 1+y
        for i,z in enumerate(zt):
            if z < 0:
                zt[i] = 1+z
        obj = Run_chain_HOP(self.period, self.padding,
            self.num_neighbors, self.bounds,
            xt,
            yt,
            zt,
            self.particle_fields["ParticleMassMsun"]/self.total_mass,
            self.threshold)
        self.densities, self.tags = obj.dens, obj.gIDs
        self.particle_fields["densities"] = self.densities
        self.particle_fields["tags"] = self.tags
        if self._distributed:
            self.padded_particles = obj.padded_particles
            self.particle_fields["particle_index"] = self._data_source["particle_index"][self._base_indices]
            mine, halo_info = self._mpi_info_dict(max(self.tags) + 1)
            nhalos = sum(halo_info.values())
            # Figure out our offset
            my_first_id = sum([v for k,v in halo_info.items() if k < mine])
            after = my_first_id + max(self.tags) + 1
            # build the list of [particle_index, threadID, grpID]
            if len(self.padded_particles):
                self.particle_index_map = self.particle_fields["particle_index"][self.padded_particles]
                tag_map = self.tags[self.padded_particles]
                temp_map = []
                for i,particle in enumerate(self.particle_index_map):
                    temp_map.append([particle,mine,tag_map[i]])
                self.particle_index_map = na.array(temp_map)
            else:
                self.particle_index_map = na.array([])
            # concatenate the map on all threads, it would be faster to
            # intelligently send this to only the processors that need each 
            # part of the data, but that's harder. Maybe later.
            #mylog.info("%s" % str(self.particle_index_map))
            self.particle_index_map = self._mpi_catarray(self.particle_index_map)
            if len(self.particle_index_map) == 0:
                # none of the groups cross the boundaries, so we're done here.
                mylog.info("No groups need to be merged.")
                return
            # but if there are crossings, move forward.
            # On each processor, if one of the particles it owns (in the real
            # region!) is in the map, add an entry to the thread_group_map
            self.thread_group_map = []
            particle_map_column = self.particle_index_map[:,0]
            for index,ownPart in enumerate(self.particle_fields["particle_index"]):
                if ownPart in particle_map_column:
                    remote_index = 0
                    while particle_map_column[remote_index] != ownPart:
                        remote_index += 1
                    point = [xt[index],
                        yt[index],
                        zt[index]]
                    if self._is_inside(point) is True:
                        #mylog.info('point %s bound %s ownPart %s' % \
                        #(str(point),str(self.bounds),str(ownPart)))
                        # find out which local group this particle lives in
                        local_groupID = self.particle_fields["tags"][index]
                        # add a mapping, if it isn't already there
                        # remote->local, remote has priority!
                        map = [self.particle_index_map[remote_index][1],\
                        self.particle_index_map[remote_index][2], \
                        mine, local_groupID]
                        if map not in self.thread_group_map:
                            self.thread_group_map.append(map)
            # concatenate the thread_group_map on all threads
            self.thread_group_map = self._mpi_catlist(self.thread_group_map)
            print self.thread_group_map

    def write_out(self, filename="chainHopAnalysis.out"):
        HaloList.write_out(self, filename)

class GenericHaloFinder(ParallelAnalysisInterface):
    def __init__(self, pf, dm_only=True, padding=0.0):
        self.pf = pf
        self.hierarchy = pf.h
        self.center = (pf["DomainRightEdge"] + pf["DomainLeftEdge"])/2.0

    def _parse_halolist(self, threshold_adjustment):
        groups, max_dens, hi  = [], {}, 0
        LE, RE = self.bounds
        for halo in self._groups:
            this_max_dens = halo.maximum_density_location()
            # if the most dense particle is in the box, keep it
            if na.all((this_max_dens >= LE) & (this_max_dens <= RE)):
                # Now we add the halo information to OURSELVES, taken from the
                # self.hop_list
                # We need to mock up the HOPHaloList thingie, so we need to set:
                #     self._max_dens
                max_dens_temp = list(self._max_dens[halo.id])[0] / threshold_adjustment
                max_dens[hi] = [max_dens_temp] + list(self._max_dens[halo.id])[1:4]
                groups.append(self._halo_class(self, hi))
                groups[-1].indices = halo.indices
                self._claim_object(groups[-1])
                hi += 1
        del self._groups, self._max_dens # explicit >> implicit
        self._groups = groups
        self._max_dens = max_dens

    def _join_halolists(self):
        # First we get the total number of halos the entire collection
        # has identified
        # Note I have added a new method here to help us get information
        # about processors and ownership and so forth.
        # _mpi_info_dict returns a dict of {proc: whatever} where whatever is
        # what is fed in on each proc.
        mine, halo_info = self._mpi_info_dict(len(self))
        nhalos = sum(halo_info.values())
        # Figure out our offset
        my_first_id = sum([v for k,v in halo_info.items() if k < mine])
        # Fix our max_dens
        max_dens = {}
        for i,m in self._max_dens.items(): max_dens[i+my_first_id] = m
        self._max_dens = max_dens
        # sort the list by the size of the groups
        # Now we add ghost halos and reassign all the IDs
        # Note: we already know which halos we own!
        after = my_first_id + len(self._groups)
        # One single fake halo, not owned, does the trick
        self._groups = [self._halo_class(self, i) for i in range(my_first_id)] + \
                       self._groups + \
                       [self._halo_class(self, i) for i in range(after, nhalos)]
        id = 0
        for proc in sorted(halo_info.keys()):
            for halo in self._groups[id:id+halo_info[proc]]:
                halo.id = id
                halo._distributed = self._distributed
                halo._owner = proc
                id += 1
        def haloCmp(h1,h2):
            c = cmp(h1.get_size(),h2.get_size())
            if c != 0:
                return -1 * c
            if c == 0:
                return cmp(h1.center_of_mass()[0],h2.center_of_mass()[0])
        self._groups.sort(haloCmp)
        sorted_max_dens = {}
        for i, halo in enumerate(self._groups):
            if halo.id in self._max_dens:
                sorted_max_dens[i] = self._max_dens[halo.id]
            halo.id = i
        self._max_dens = sorted_max_dens
        
    def _reposition_particles(self, bounds):
        # This only does periodicity.  We do NOT want to deal with anything
        # else.  The only reason we even do periodicity is the 
        LE, RE = bounds
        dw = self.pf["DomainRightEdge"] - self.pf["DomainLeftEdge"]
        for i, ax in enumerate('xyz'):
            arr = self._data_source["particle_position_%s" % ax]
            arr[arr < LE[i]-self.padding] += dw[i]
            arr[arr > RE[i]+self.padding] -= dw[i]

    def write_out(self, filename):
        self._data_source.get_data(["particle_velocity_%s" % ax for ax in 'xyz'])
        f = self._write_on_root(filename)
        HaloList.write_out(self, f)

    @parallel_blocking_call
    def write_particle_lists(self, prefix):
        fn = "%s.h5" % self._get_filename(prefix)
        f = h5py.File(fn, "w")
        for halo in self._groups:
            if not self._is_mine(halo): continue
            halo.write_particle_list(f)

class chainHF(GenericHaloFinder, chainHOPHaloList):
    def __init__(self, pf, threshold=160, dm_only=True):
        GenericHaloFinder.__init__(self, pf, dm_only, padding=0.0)
        self.padding = 0.0 
        self.num_neighbors = 65
        # get the total number of particles across all procs, with no padding
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        # also get the total mass of particles
        if dm_only:
            select = self._data_source["creation_time"] < 0
            total_mass = self._mpi_allsum((self._data_source["ParticleMassMsun"][select]).sum())
            n_parts = self._mpi_allsum((self._data_source["particle_position_x"][select]).size)
        else:
            total_mass = self._mpi_allsum(self._data_source["ParticleMassMsun"].sum())
            n_parts = self._mpi_allsum(self._data_source["particle_position_x"].size)
        # get the average spacing between particles
        l = pf["DomainRightEdge"] - pf["DomainLeftEdge"]
        vol = l[0] * l[1] * l[2]
        avg_spacing = (float(vol) / n_parts)**(1./3.)
        # padding is a function of inter-particle spacing, this is an
        # approximation, but it's OK with the safety factor
        safety = 4
        self.padding = (self.num_neighbors)**(1./3.) * safety * avg_spacing
        print 'padding',self.padding,'avg_spacing',avg_spacing,'vol',vol,'nparts',n_parts
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        if not padded:
            self.padding = 0.0
        self.bounds = (LE, RE)
        chainHOPHaloList.__init__(self, self._data_source, self.padding, \
        self.num_neighbors, self.bounds, total_mass, threshold=threshold, dm_only=dm_only)
        self._parse_halolist(1.)
        self._join_halolists()

class HOPHaloFinder(GenericHaloFinder, HOPHaloList):
    def __init__(self, pf, threshold=160, dm_only=True, padding=0.02):
        GenericHaloFinder.__init__(self, pf, dm_only, padding)
        
        # do it once with no padding so the total_mass is correct (no duplicated particles)
        self.padding = 0.0
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        # For scaling the threshold, note that it's a passthrough
        if dm_only:
            select = self._data_source["creation_time"] > 0
            total_mass = self._mpi_allsum((self._data_source["ParticleMassMsun"][select]).sum())
            sub_mass = (self._data_source["ParticleMassMsun"][select]).sum()
        else:
            total_mass = self._mpi_allsum(self._data_source["ParticleMassMsun"].sum())
            sub_mass = self._data_source["ParticleMassMsun"].sum()
        # MJT: Note that instead of this, if we are assuming that the particles
        # are all on different processors, we should instead construct an
        # object representing the entire domain and sum it "lazily" with
        # Derived Quantities.
        self.padding = padding #* pf["unitary"] # This should be clevererer
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        self.bounds = (LE, RE)
        # reflect particles around the periodic boundary
        #self._reposition_particles((LE, RE))
        #sub_mass = self._data_source["ParticleMassMsun"].sum()
        HOPHaloList.__init__(self, self._data_source, threshold*total_mass/sub_mass, dm_only)
        self._parse_halolist(total_mass/sub_mass)
        self._join_halolists()

class FOFHaloFinder(GenericHaloFinder, FOFHaloList):
    def __init__(self, pf, link=0.2, dm_only=True, padding=0.02):
        self.pf = pf
        self.hierarchy = pf.h
        self.center = (pf["DomainRightEdge"] + pf["DomainLeftEdge"])/2.0
        self.padding = 0.0 #* pf["unitary"] # This should be clevererer
        # get the total number of particles across all procs, with no padding
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        n_parts = self._mpi_allsum(self._data_source["particle_position_x"].size)
        # get the average spacing between particles
        l = pf["DomainRightEdge"] - pf["DomainLeftEdge"]
        vol = l[0] * l[1] * l[2]
        avg_spacing = (float(vol) / n_parts)**(1./3.)
        self.padding = padding
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        self.bounds = (LE, RE)
        # reflect particles around the periodic boundary
        #self._reposition_particles((LE, RE))
        # here is where the FOF halo finder is run
        FOFHaloList.__init__(self, self._data_source, link * avg_spacing, dm_only)
        self._parse_halolist(1.)
        self._join_halolists()

HaloFinder = HOPHaloFinder

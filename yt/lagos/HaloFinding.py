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
try:
    from yt.lagos.parallelHOP.parallelHOP import *
except ImportError:
    mylog.error("ParallelHOP not imported.")

try:
    from yt.lagos.fof.EnzoFOF import RunFOF
except ImportError:
    pass
from yt.performance_counters import yt_counters, time_function

from kd import *
import math
from collections import defaultdict

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

    def __init__(self, halo_list, id, indices = None, size=None, CoM=None,
        max_dens_point=None, group_total_mass=None, max_radius=None, bulk_vel=None,
        tasks=None):
        self.halo_list = halo_list
        self.id = id
        self.data = halo_list._data_source
        if indices is not None: self.indices = halo_list._base_indices[indices]
        # We assume that if indices = None, the instantiator has OTHER plans
        # for us -- i.e., setting it somehow else
        self.size = size
        self.CoM = CoM
        self.max_dens_point = max_dens_point
        self.group_total_mass = group_total_mass
        self.max_radius = max_radius
        self.bulk_vel = bulk_vel
        self.tasks = tasks

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

class parallelHOPHalo(Halo,ParallelAnalysisInterface):
    dont_wrap = ["maximum_density","maximum_density_location",
        "center_of_mass","total_mass","bulk_velocity","maximum_radius",
        "get_size","get_sphere", "write_particle_list","__getitem__"]

    def __getitem__(self, key):
        return self.data[key][self.indices]

    def maximum_density(self):
        """
        Return the HOP-identified maximum density.
        """
        if self.max_dens_point is not None:
            return self.halo_list._max_dens[self.id][0]
        max = self._mpi_allmax(self.halo_list._max_dens[self.id][0])
        return max

    def maximum_density_location(self):
        """
        Return the location HOP identified as maximally dense.
        """
        if self.max_dens_point is not None:
            return na.array([self.halo_list._max_dens[self.id][1], self.halo_list._max_dens[self.id][2],
                self.halo_list._max_dens[self.id][3]])
        # If I own the maximum density, my location is globally correct.
        max_dens = self.maximum_density()
        if self.halo_list._max_dens[self.id][0] == max_dens:
            value = na.array([
                self.halo_list._max_dens[self.id][1],
                self.halo_list._max_dens[self.id][2],
                self.halo_list._max_dens[self.id][3]])
        else:
            value = na.array([0,0,0])
        # This works, and isn't appropriate but for now will be fine...
        value = self._mpi_allsum(value)
        return value

    def center_of_mass(self):
        """
        Calculate and return the center of mass.
        """
        # If it's precomputed, we save time!
        if self.CoM is not None:
            return self.CoM
        # This need to be called by all tasks, but not all will end up using
        # it.
        c_vec = self.maximum_density_location() - na.array([0.5,0.5,0.5])
        if self.indices is not None:
            pm = self["ParticleMassMsun"]
            cx = (self["particle_position_x"] - c_vec[0])
            cy = (self["particle_position_y"] - c_vec[1])
            cz = (self["particle_position_z"] - c_vec[2])
            com = na.array([v-na.floor(v) for v in [cx,cy,cz]])
            my_mass = pm.sum()
            my_com = ((com*pm).sum(axis=1)/my_mass + c_vec) * my_mass
        else:
            my_mass = 0.
            my_com = na.array([0.,0.,0.])
        global_mass = self._mpi_allsum(my_mass)
        global_com = self._mpi_allsum(my_com)
        return global_com / global_mass

    def total_mass(self):
        """
        Returns the total mass in solar masses of the halo.
        """
        if self.group_total_mass is not None:
            return self.group_total_mass
        if self.indices is not None:
            my_mass = self["ParticleMassMsun"].sum()
        else:
            my_mass = 0.
        global_mass = self._mpi_allsum(float(my_mass))
        return global_mass

    def bulk_velocity(self):
        """
        Returns the mass-weighted average velocity.
        """
        if self.bulk_vel is not None:
            return self.bulk_vel
        # Unf. this cannot be reasonably computed inside of parallelHOP because
        # we don't pass velocities in.
        if self.indices is not None:
            pm = self["ParticleMassMsun"]
            vx = (self["particle_velocity_x"] * pm).sum()
            vy = (self["particle_velocity_y"] * pm).sum()
            vz = (self["particle_velocity_z"] * pm).sum()
            pm = pm.sum()
        else:
            pm = 0.
            vx = 0.
            vy = 0.
            vz = 0.
        bv = na.array([vx,vy,vz,pm])
        global_bv = self._mpi_allsum(bv)
        return global_bv[:3]/global_bv[3]

    def maximum_radius(self, center_of_mass=True):
        """
        Returns the maximum radius in the halo for all particles,
        either from the point of maximum density or from the (default)
        *center_of_mass*.
        """
        if self.max_radius is not None:
            return self.max_radius
        if center_of_mass: center = self.center_of_mass()
        else: center = self.maximum_density_location()
        if self.indices is not None:
            rx = na.abs(self["particle_position_x"]-center[0])
            ry = na.abs(self["particle_position_y"]-center[1])
            rz = na.abs(self["particle_position_z"]-center[2])
            r = na.sqrt(na.minimum(rx, 1.0-rx)**2.0
                    +   na.minimum(ry, 1.0-ry)**2.0
                    +   na.minimum(rz, 1.0-rz)**2.0)
            my_max = r.max()
            
        else:
            my_max = 0.
        return self._mpi_allmax(my_max)

    def get_size(self):
        if self.size is not None:
            return self.size
        if self.indices is not None:
            my_size = self.indices.size
        else:
            my_size = 0
        global_size = self._mpi_allsum(my_size)
        return global_size


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

    def write_particle_lists_txt(self, prefix, fp=None):
        """
        Write out the location of halo data in hdf5 files to *prefix*.
        """
        if hasattr(fp, 'write'):
            f = fp
        else:
            f = open("%s.txt" % prefix,"w")
        for group in self:
            if group.tasks is not None:
                fn = ""
                for task in group.tasks:
                    fn += "%s.h5 " % self._get_filename(prefix, rank=task)
            elif self._distributed:
                fn = "%s.h5" % self._get_filename(prefix, rank=group._owner)
            else:
                fn = "%s.h5" % self._get_filename(prefix)
            gn = "Halo%08i" % (group.id)
            f.write("%s %s\n" % (gn, fn))
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

class parallelHOPHaloList(HaloList,ParallelAnalysisInterface):
    _name = "parallelHOP"
    _halo_class = parallelHOPHalo
    _fields = ["particle_position_%s" % ax for ax in 'xyz'] + \
              ["ParticleMassMsun", "particle_index"]
    
    def __init__(self, data_source, padding, num_neighbors, bounds, total_mass,
        period, threshold=160.0, dm_only=True, rearrange=True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self.threshold = threshold
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.total_mass = total_mass
        self.rearrange = rearrange
        self.period = period
        mylog.info("Initializing HOP")
        HaloList.__init__(self, data_source, dm_only)

    def _run_finder(self):
        yt_counters("Reading Data")
        obj = RunParallelHOP(self.period, self.padding,
            self.num_neighbors, self.bounds,
            self.particle_fields["particle_position_x"],
            self.particle_fields["particle_position_y"],
            self.particle_fields["particle_position_z"],
            self.particle_fields["particle_index"],
            self.particle_fields["ParticleMassMsun"]/self.total_mass,
            self.threshold, rearrange=self.rearrange)
        self.densities, self.tags = obj.density, obj.chainID
        self.group_count = obj.group_count
        self.group_sizes = obj.group_sizes
        self.CoM = obj.CoM
        self.Tot_M = obj.Tot_M * self.total_mass
        self.max_dens_point = obj.max_dens_point
        self.max_radius = obj.max_radius
        # Precompute the bulk velocity in parallel.
        yt_counters("Precomp bulk vel.")
        self.bulk_vel = na.zeros((self.group_count, 3), dtype='float64')
        yt_counters("bulk vel. reading data")
        pm = self._data_source["ParticleMassMsun"][self._base_indices]
        xv = self._data_source["particle_velocity_x"][self._base_indices]
        yv = self._data_source["particle_velocity_y"][self._base_indices]
        zv = self._data_source["particle_velocity_z"][self._base_indices]
        yt_counters("bulk vel. reading data")
        yt_counters("bulk vel. computing")
        for i, pmass in enumerate(pm):
            groupID = self.tags[i]
            if groupID == -1: continue
            self.bulk_vel[groupID] += (na.array([xv[i], yv[i], zv[i]])*pmass)
        # Bring it together, and divide by the previously computed total mass
        # of each halo.
        self.bulk_vel = self._mpi_Allsum_double(self.bulk_vel)
        for groupID in xrange(self.group_count):
            self.bulk_vel[groupID] = self.bulk_vel[groupID] / self.Tot_M[groupID]
        yt_counters("bulk vel. computing")
        self.taskID = obj.mine
        self.halo_taskmap = obj.halo_taskmap # A defaultdict.
        del obj
        yt_counters("Precomp bulk vel.")

    def _parse_output(self):
        yt_counters("Final Grouping")
        """
        Each task will make an entry for all groups, but it may be empty.
        """
        unique_ids = na.unique(self.tags)
        counts = na.bincount((self.tags+1).tolist())
        sort_indices = na.argsort(self.tags)
        grab_indices = na.indices(self.tags.shape).ravel()[sort_indices]
        dens = self.densities[sort_indices]
        cp = 0
        index = 0
        for i in unique_ids:
            if i == -1:
                cp += counts[i+1]
                continue
            # If there is a gap in the unique_ids, make empty groups to 
            # fill it in.
            while index < i:
                self._groups.append(self._halo_class(self, index, \
                    size=self.group_sizes[index], CoM=self.CoM[index], \
                    max_dens_point=self.max_dens_point[index], \
                    group_total_mass=self.Tot_M[index], max_radius=self.max_radius[index],
                    bulk_vel=self.bulk_vel[index], tasks=self.halo_taskmap[index]))
                # I don't own this halo
                self._do_not_claim_object(self._groups[-1])
                self._max_dens[index] = (self.max_dens_point[index][0], self.max_dens_point[index][1], \
                    self.max_dens_point[index][2], self.max_dens_point[index][3])
                index += 1
            cp_c = cp + counts[i+1]
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(self._halo_class(self, i, group_indices, \
                size=self.group_sizes[i], CoM=self.CoM[i], \
                max_dens_point=self.max_dens_point[i], \
                group_total_mass=self.Tot_M[i], max_radius=self.max_radius[i],
                bulk_vel=self.bulk_vel[i], tasks=self.halo_taskmap[index]))
            # This halo may be owned by many, including this task
            self._claim_object(self._groups[-1])
            self._max_dens[i] = (self.max_dens_point[i][0], self.max_dens_point[i][1], \
                self.max_dens_point[i][2], self.max_dens_point[i][3])
            cp += counts[i+1]
            index += 1
        # If there are missing groups at the end, add them.
        while index < self.group_count:
            self._groups.append(self._halo_class(self, index, \
            size=self.group_sizes[index], CoM=self.CoM[index], \
            max_dens_point=self.max_dens_point[i], \
            group_total_mass=self.Tot_M[index], max_radius=self.max_radius[index],
            bulk_vel=self.bulk_vel[index], tasks=self.halo_taskmap[index]))
            self._do_not_claim_object(self._groups[-1])
            self._max_dens[index] = (self.max_dens_point[index][0], self.max_dens_point[index][1], \
                self.max_dens_point[index][2], self.max_dens_point[index][3])
            index += 1

    def __len__(self):
        return self.group_count

    def write_out(self, filename="parallelHopAnalysis.out"):
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

    def write_particle_lists_txt(self, prefix):
        f = self._write_on_root("%s.txt" % prefix)
        HaloList.write_particle_lists_txt(self, prefix, fp=f)

    @parallel_blocking_call
    def write_particle_lists(self, prefix):
        fn = "%s.h5" % self._get_filename(prefix)
        f = h5py.File(fn, "w")
        for halo in self._groups:
            if not self._is_mine(halo): continue
            halo.write_particle_list(f)

class parallelHF(GenericHaloFinder, parallelHOPHaloList):
    def __init__(self, pf, threshold=160, dm_only=True, resize=False, rearrange=True,\
        fancy_padding=True, safety=2.5):
        GenericHaloFinder.__init__(self, pf, dm_only, padding=0.0)
        self.padding = 0.0 
        self.num_neighbors = 65
        self.safety = safety
        period = pf["DomainRightEdge"] - pf["DomainLeftEdge"]
        # get the total number of particles across all procs, with no padding
        padded, LE, RE, self._data_source = self._partition_hierarchy_3d(padding=self.padding)
        # also get the total mass of particles
        yt_counters("Reading Data")
        # Adaptive subregions by bisection.
        ds_names = ["particle_position_x","particle_position_y","particle_position_z"]
        if resize and self._mpi_get_size()!=None:
            cut_list = self._partition_hierarchy_3d_bisection_list()
            for i,cut in enumerate(cut_list):
                dim = cut[0]
                if i == 0:
                    width = pf["DomainRightEdge"][dim] - pf["DomainLeftEdge"][dim]
                    new_LE = pf["DomainLeftEdge"]
                else:
                    new_LE, new_RE = new_top_bounds
                    width = new_RE[dim] - new_LE[dim]
                data = self._data_source[ds_names[dim]]
                if i == 0:
                    local_parts = data.size
                    n_parts = self._mpi_allsum(local_parts)
                    min = self._mpi_allmin(local_parts)
                    max = self._mpi_allmax(local_parts)
                    mylog.info("Initial distribution: (min,max,tot) (%d,%d,%d)" \
                        % (min, max, n_parts))
                num_bins = 1000
                bin_width = float(width)/float(num_bins)
                bins = na.arange(num_bins+1, dtype='float64') * bin_width + new_LE[dim]
                counts, bins = na.histogram(data, bins, new=True)
                if i == 0:
                    new_group, new_comm, LE, RE, new_top_bounds, new_cc, self._data_source= \
                        self._partition_hierarchy_3d_bisection(dim, bins, counts, top_bounds = None,\
                        old_group = None, old_comm = None, cut=cut)
                else:
                    new_group, new_comm, LE, RE, new_top_bounds, new_cc, self._data_source = \
                        self._partition_hierarchy_3d_bisection(dim, bins, counts, top_bounds = new_top_bounds,\
                        old_group = new_group, old_comm = new_comm, cut=cut, old_cc=new_cc)
        # get the average spacing between particles for this region
        # The except is for the serial case, where the full box is what we want.
        data = self._data_source["particle_position_x"]
        try:
            l = self._data_source.right_edge - self._data_source.left_edge
        except AttributeError:
            l = pf["DomainRightEdge"] - pf["DomainLeftEdge"]
        vol = l[0] * l[1] * l[2]
        full_vol = vol
        if not fancy_padding:
            avg_spacing = (float(vol) / data.size)**(1./3.)
            # padding is a function of inter-particle spacing, this is an
            # approximation, but it's OK with the safety factor
            padding = (self.num_neighbors)**(1./3.) * self.safety * avg_spacing
            self.padding = (na.ones(3,dtype='float64')*padding, na.ones(3,dtype='float64')*padding)
            mylog.info('padding %s avg_spacing %f vol %f local_parts %d' % \
                (str(self.padding), avg_spacing, vol, data.size))
        # Another approach to padding, perhaps more accurate.
        elif fancy_padding and self._distributed:
            LE_padding, RE_padding = na.empty(3,dtype='float64'), na.empty(3,dtype='float64')
            for dim in xrange(3):
                data = self._data_source[ds_names[dim]]
                num_bins = 1000
                width = self._data_source.right_edge[dim] - self._data_source.left_edge[dim]
                area = (self._data_source.right_edge[(dim+1)%3] - self._data_source.left_edge[(dim+1)%3]) * \
                    (self._data_source.right_edge[(dim+2)%3] - self._data_source.left_edge[(dim+2)%3])
                bin_width = float(width)/float(num_bins)
                bins = na.arange(num_bins+1, dtype='float64') * bin_width + self._data_source.left_edge[dim]
                counts, bins = na.histogram(data, bins, new=True)
                # left side.
                start = 0
                count = counts[0]
                while count < self.num_neighbors:
                    start += 1
                    count += counts[start]
                # Get the avg spacing in just this boundary.
                vol = area * (bins[start+1] - bins[0])
                avg_spacing = (float(vol) / count)**(1./3.)
                LE_padding[dim] = (self.num_neighbors)**(1./3.) * self.safety * avg_spacing
                # right side.
                start = -1
                count = counts[-1]
                while count < self.num_neighbors:
                    start -= 1
                    count += counts[start]
                vol = area * (bins[-1] - bins[start-1])
                avg_spacing = (float(vol) / count)**(1./3.)
                RE_padding[dim] = (self.num_neighbors)**(1./3.) * self.safety * avg_spacing
            self.padding = (LE_padding, RE_padding)
            mylog.info('fancy_padding %s avg_spacing %f full_vol %f local_parts %d %s' % \
                (str(self.padding), avg_spacing, full_vol, data.size, str(self._data_source)))
        # Now we get the full box mass after we have the final composition of
        # subvolumes.
        total_mass = self._mpi_allsum(self._data_source["ParticleMassMsun"].sum())
        if not self._distributed:
            self.padding = (na.zeros(3,dtype='float64'), na.zeros(3,dtype='float64'))
        self.bounds = (LE, RE)
        (LE_padding, RE_padding) = self.padding
        parallelHOPHaloList.__init__(self, self._data_source, self.padding, \
        self.num_neighbors, self.bounds, total_mass, period, threshold=threshold, dm_only=dm_only, rearrange=rearrange)
        self._join_halolists()
        yt_counters("Final Grouping")

    def _join_halolists(self):
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
